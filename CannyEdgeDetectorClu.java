/**
 * This program performs the Canny Edge Detection in a cluster parallel computer
 * It mainly detects the prominent edges in the image and display then in the output.
 * The images are sliced by rows and send to multiple computers and then reduced to final result.
 * 
 * @author Divin Visariya
 * @author Sreeprasad Govindankutty
 * @author Varun Goyal
 * 
 */

import java.awt.Image;
import java.awt.image.BufferedImage;
import java.awt.image.ImageObserver;
import java.awt.image.PixelGrabber;
import java.io.File;
import java.io.InputStream;
import java.util.Arrays;

import javax.imageio.ImageIO;
import edu.rit.image.PJGColorImage;
import edu.rit.mp.IntegerBuf;
import edu.rit.pj.Comm;
import edu.rit.util.Range;

import edu.rit.pj.WorkerIntegerForLoop;
import edu.rit.pj.WorkerRegion;
import edu.rit.pj.WorkerTeam;
import edu.rit.pj.reduction.IntegerOp;
 

public class CannyEdgeDetectorClu {

	static Comm world;
	static int size;
	static int rank;
	static InputStream is;
	static PJGColorImage old_image;
	static PJGColorImage new_image;

	static int initX, maxX, initY, maxY;

	static IntegerBuf[] slices;
	static IntegerBuf myslice;

	// Image matrix.
	static int[][] matrix;
	static BufferedImage image;
	static Range[] ranges;
	static Range myrange;
	static int mylb;
	static int myub;
	// statics
	private final static float GAUSSIAN_CUT_OFF = 0.005f;
	private final static float MAGNITUDE_SCALE = 100F;
	private final static float MAGNITUDE_LIMIT = 1000F;
	private final static int MAGNITUDE_MAX = (int) (MAGNITUDE_SCALE * MAGNITUDE_LIMIT);

	// fields
	static private int pixels[];
	private static int height;
	private static int width;
	private static int picsize;
	private static int[] data;
	private static int[] magnitude;
	private static BufferedImage sourceImage;
	private static BufferedImage edgesImage;

	private static float gaussianKernelRadius = 2f;
	private static float lowThreshold = 2.5f;
	private static float highThreshold = 7.5f;
	private static int gaussianKernelWidth = 16;
	private static boolean contrastNormalized = true;

	private static float[] xConv;
	private static float[] yConv;
	private static float[] xGradient;
	private static float[] yGradient;

	
	private static long start_time,end_time;


	// accessors

	/**
	 * Specifies the image that will provide the luminance data in which edges
	 * will be detected. A source image must be set before the process method is
	 * called.
	 * 
	 * @param image    a source of luminance data
	 */

	public static void setSourceImage(BufferedImage image) {
		sourceImage = image;
	}

	/**
	 * Obtains an image containing the edges detected during the last call to
	 * the process method. The buffered image is an opaque image of type
	 * BufferedImage.TYPE_INT_ARGB in which edge pixels are white and all other
	 * pixels are black.
	 * 
	 * @return an image containing the detected edges, or null if the process
	 *         method has not yet been called.
	 */

	public static BufferedImage getEdgesImage() {
		return edgesImage;
	}

	/**
	 * Sets the low threshold for hysteresis. Suitable values for this parameter
	 * must be determined experimentally for each application. It is nonsensical
	 * (though not prohibited) for this value to exceed the high threshold
	 * value.
	 * 
	 * @param threshold
	 *            a low hysteresis threshold
	 */

	public static void setLowThreshold(float threshold) {
		if (threshold < 0)
			throw new IllegalArgumentException();
		lowThreshold = threshold;
	}

	/**
	 * Sets the high threshold for hysteresis. Suitable values for this
	 * parameter must be determined experimentally for each application. It is
	 * nonsensical (though not prohibited) for this value to be less than the
	 * low threshold value.
	 * 
	 * @param threshold
	 *            a high hysteresis threshold
	 */

	public static void setHighThreshold(float threshold) {
		if (threshold < 0)
			throw new IllegalArgumentException();
		highThreshold = threshold;
	}

	/**
	 * The number of pixels across which the Gaussian kernel is applied. This
	 * implementation will reduce the radius if the contribution of pixel values
	 * is deemed negligable, so this is actually a maximum radius.
	 * 
	 * @param gaussianKernelWidth
	 *            a radius for the convolution operation in pixels, at least 2.
	 */

	public void setGaussianKernelWidth(int gaussianKernelWidth) {
		if (gaussianKernelWidth < 2)
			throw new IllegalArgumentException();
		CannyEdgeDetectorClu.gaussianKernelWidth = gaussianKernelWidth;
	}

	/**
	 * Sets the radius of the Gaussian convolution kernel used to smooth the
	 * source image prior to gradient calculation.
	 * 
	 * @return a Gaussian kernel radius in pixels, must exceed 0.1f.
	 */

	public void setGaussianKernelRadius(float gaussianKernelRadius) {
		if (gaussianKernelRadius < 0.1f)
			throw new IllegalArgumentException();
		CannyEdgeDetectorClu.gaussianKernelRadius = gaussianKernelRadius;
	}

	/**
	 * Whether the luminance data extracted from the source image is normalized
	 * by linearizing its histogram prior to edge extraction. The default value
	 * is false.
	 * 
	 * @return whether the contrast is normalized
	 */

	public boolean isContrastNormalized() {
		return contrastNormalized;
	}

	/**
	 * Sets whether the contrast is normalized
	 * 
	 * @param contrastNormalized
	 *            true if the contrast should be normalized, false otherwise
	 */

	public void setContrastNormalized(boolean contrastNormalized) {
		CannyEdgeDetectorClu.contrastNormalized = contrastNormalized;
	}

	
	/**
	 * The main processing is done by this method which calls all other steps of Canny edge detection
	 * @param w (width of the image)
	 * @throws Exception
	 */
	public static void process(int w) throws Exception {
		width = w;
		picsize = width * height;
		initArrays();
		readLuminance();
		if (contrastNormalized)
			normalizeContrast();
		computeGradients();
		int low = Math.round(lowThreshold * MAGNITUDE_SCALE);
		int high = Math.round(highThreshold * MAGNITUDE_SCALE);
		performHysteresis(low, high);
		thresholdEdges();

		// All process results will be send to process 0
		// for reduction
		IntegerBuf buf = IntegerBuf.buffer(data);
		world.reduce(0, buf, IntegerOp.SUM);

	}

	
	/**
	 * Initializes the arrays
	 */
	private static void initArrays() {
		if (data == null || picsize != data.length) {
			data = new int[picsize];
			magnitude = new int[picsize];

			xConv = new float[picsize];
			yConv = new float[picsize];
			xGradient = new float[picsize];
			yGradient = new float[picsize];
		}
	}
	
	
	/**
	 * This computes the gradient values for the image 
	 * It mainly uses the gaussian mask to compute the values
	 * 
	 * @throws Exception
	 */
	private static void computeGradients()
			throws Exception {

		int kernelWidth= gaussianKernelWidth;
		float kernelRadius=gaussianKernelRadius;
		
		float kernel[] = new float[kernelWidth];
		float diffKernel[] = new float[kernelWidth];
		int kwidth;
		for (kwidth = 0; kwidth < kernelWidth; kwidth++) {
			float g1 = gaussian(kwidth, kernelRadius);
			if (g1 <= GAUSSIAN_CUT_OFF && kwidth >= 2)
				break;
			float g2 = gaussian(kwidth - 0.5f, kernelRadius);
			float g3 = gaussian(kwidth + 0.5f, kernelRadius);
			kernel[kwidth] = (g1 + g2 + g3) / 3f
					/ (2f * (float) Math.PI * kernelRadius * kernelRadius);
			diffKernel[kwidth] = g3 - g2;
		}

		int initX = kwidth - 1;
		int maxX = width - (kwidth - 1);
		int initY = width * (kwidth - 1);
		int maxY = width * (height - (kwidth - 1));

		// perform convolution in x and y directions
		for (int x = initX; x < maxX; x++) {
			for (int y = initY; y < maxY; y += width) {
				int index = x + y;
				float sumX = data[index] * kernel[0];
				float sumY = sumX;
				int xOffset = 1;
				int yOffset = width;
				for (; xOffset < kwidth;) {
					sumY += kernel[xOffset]
							* (data[index - yOffset] + data[index + yOffset]);
					sumX += kernel[xOffset]
							* (data[index - xOffset] + data[index + xOffset]);
					yOffset += width;
					xOffset++;
				}

				yConv[index] = sumY;
				xConv[index] = sumX;
			}

		}

		//xGradient is computed here
		for (int x = initX; x < maxX; x++) {
			for (int y = initY; y < maxY; y += width) {
				float sum = 0f;
				int index = x + y;
				for (int i = 1; i < kwidth; i++)
					sum += diffKernel[i]
							* (yConv[index - i] - yConv[index + i]);

				xGradient[index] = sum;
			}

		}

		//yGradient is computed here
		for (int x = kwidth; x < width - kwidth; x++) {
			for (int y = initY; y < maxY; y += width) {
				float sum = 0.0f;
				int index = x + y;
				int yOffset = width;
				for (int i = 1; i < kwidth; i++) {
					sum += diffKernel[i]
							* (xConv[index - yOffset] - xConv[index + yOffset]);
					yOffset += width;
				}

				yGradient[index] = sum;
			}

		}

		initX = kwidth;
		maxX = width - kwidth;
		initY = width * kwidth;
		maxY = width * (height - kwidth);

		final int init_x = initX;
		final int init_y = initY;
		final int max_x = maxX;
		final int max_y = maxY;
		
		//This does the computation for all 8 directions from the pixel
		new WorkerTeam().execute(new WorkerRegion() {
			public void run() throws Exception {
				execute(init_x, max_x - 2, new WorkerIntegerForLoop() {
					public void run(int first, int last) {
						for (int x = first; x < last; x++) {
							for (int y = init_y; y < max_y; y += width) {
								
								//Variables storing the eigth pixels surrounding the pixel under consideration
								int index = x + y;
								int indexN = index - width;
								int indexS = index + width;
								int indexW = index - 1;
								int indexE = index + 1;
								int indexNW = indexN - 1;
								int indexNE = indexN + 1;
								int indexSW = indexS - 1;
								int indexSE = indexS + 1;

								float xGrad = xGradient[index];
								float yGrad = yGradient[index];
								float gradMag = hypot(xGrad, yGrad);

								// perform non-maximal supression
								// It calculates the sqrt of sum of xgradient^2 and ygradient^2
								float nMag = hypot(xGradient[indexN],
										yGradient[indexN]);
								float sMag = hypot(xGradient[indexS],
										yGradient[indexS]);
								float wMag = hypot(xGradient[indexW],
										yGradient[indexW]);
								float eMag = hypot(xGradient[indexE],
										yGradient[indexE]);
								float neMag = hypot(xGradient[indexNE],
										yGradient[indexNE]);
								float seMag = hypot(xGradient[indexSE],
										yGradient[indexSE]);
								float swMag = hypot(xGradient[indexSW],
										yGradient[indexSW]);
								float nwMag = hypot(xGradient[indexNW],
										yGradient[indexNW]);
								float tmp;

								/*
								 * This performs the "non-maximal supression" phase of 
								 * the Canny Edge Detection in which we
								 * need to compare the gradient magnitude to that in the
								 * direction of the gradient; only if the value is a local
								 * maximum do we consider the point as an edge candidate.
								 * 
								 * We need to break the comparison into a number of different
								 * cases depending on the gradient direction so that the
								 * appropriate values can be used. To avoid computing the
								 * gradient direction, we use two simple comparisons: first we
								 * check that the partial derivatives have the same sign (1)
								 * and then we check which is larger (2). As a consequence, we
								 * have reduced the problem to one of four identical cases that
								 * each test the central gradient magnitude against the values at
								 * two points with 'identical support'; what this means is that
								 * the geometry required to accurately interpolate the magnitude
								 * of gradient function at those points has an identical
								 * geometry (upto right-angled-rotation/reflection).
								 * 
								 * When comparing the central gradient to the two interpolated
								 * values, we avoid performing any divisions by multiplying both
								 * sides of each inequality by the greater of the two partial
								 * derivatives. The common comparison is stored in a temporary
								 * variable (3) and reused in the mirror case (4).
								 * 
								 */
								if (xGrad * yGrad <= (float) 0 /* (1) */
								? Math.abs(xGrad) >= Math.abs(yGrad) /* (2) */
								? (tmp = Math.abs(xGrad * gradMag)) >= Math
										.abs(yGrad * neMag - (xGrad + yGrad)
												* eMag) /* (3) */
										&& tmp > Math.abs(yGrad * swMag
												- (xGrad + yGrad) * wMag) /* (4) */
								: (tmp = Math.abs(yGrad * gradMag)) >= Math
										.abs(xGrad * neMag - (yGrad + xGrad)
												* nMag) /* (3) */
										&& tmp > Math.abs(xGrad * swMag
												- (yGrad + xGrad) * sMag) /* (4) */
								: Math.abs(xGrad) >= Math.abs(yGrad) /* (2) */
								? (tmp = Math.abs(xGrad * gradMag)) >= Math
										.abs(yGrad * seMag + (xGrad - yGrad)
												* eMag) /* (3) */
										&& tmp > Math.abs(yGrad * nwMag
												+ (xGrad - yGrad) * wMag) /* (4) */
								: (tmp = Math.abs(yGrad * gradMag)) >= Math
										.abs(xGrad * seMag + (yGrad - xGrad)
												* sMag) /* (3) */
										&& tmp > Math.abs(xGrad * nwMag
												+ (yGrad - xGrad) * nMag) /* (4) */
								) {
									magnitude[index] = gradMag >= MAGNITUDE_LIMIT ? MAGNITUDE_MAX
											: (int) (MAGNITUDE_SCALE * gradMag);
								} else {
									magnitude[index] = 0;
								}
							}
						}
					}
				});

			}
		});
	}

	/**
	 * The hypotenunse is calculated here
	 * @param x
	 * @param y
	 * @return
	 */
	private static float hypot(float x, float y) {
		return (float) Math.hypot(x, y);
	}

	//Gaussian formula 
	private static float gaussian(float x, float sigma) {
		return (float) Math.exp(-(x * x) / (2f * sigma * sigma));
	}

	/**
	 * We check the values to get the prominent edges between low and high threshold
	 * 
	 * @param low
	 * @param high
	 */
	private static void performHysteresis(int low, int high) {

		Arrays.fill(data, 0);

		int offset = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (data[offset] == 0 && magnitude[offset] >= high) {
					follow(x, y, offset, low);
				}
				offset++;
			}
		}
	}

	/**
	 * Here we mainly follow the lines according to the threashold values
	 * All pixels with magnitude above the treshhold are considered
	 * 
	 * @param x1
	 * @param y1
	 * @param i1
	 * @param threshold
	 */
	private static void follow(int x1, int y1, int i1, int threshold) {
		int x0 = x1 == 0 ? x1 : x1 - 1;
		int x2 = x1 == width - 1 ? x1 : x1 + 1;
		int y0 = y1 == 0 ? y1 : y1 - 1;
		int y2 = y1 == height - 1 ? y1 : y1 + 1;

		data[i1] = magnitude[i1];
		for (int x = x0; x <= x2; x++) {
			for (int y = y0; y <= y2; y++) {
				int i2 = x + y * width;
				if ((y != y1 || x != x1) && data[i2] == 0
						&& magnitude[i2] >= threshold) {
					follow(x, y, i2, threshold);
					return;
				}
			}
		}
	}

	//This method converts and makes the sharp white edges on the black
	private static void thresholdEdges() {
		for (int i = 0; i < picsize; i++) {
			data[i] = data[i] > 0 ? -1 : 0xff000000;
		}
	}

	//Luminance formula for rgb value
	private static int luminance(float r, float g, float b) {
		return Math.round(0.299f * r + 0.587f * g + 0.114f * b);
	}

	/**
	 * This methods calculated the luminance of every pixel of the image.
	 * 
	 * The image are selected according to their types and luminance is 
	 * performed accordingly
	 * 
	 */
	private static void readLuminance() {
		int type = sourceImage.getType();
		
		//For RGB image
		if (type == BufferedImage.TYPE_INT_RGB
				|| type == BufferedImage.TYPE_INT_ARGB) {
			int[] pixels = (int[]) sourceImage.getData().getDataElements(0, 0,
					width, height, null);
			for (int i = 0; i < picsize; i++) {
				int p = pixels[i];
				int r = (p & 0xff0000) >> 16;
				int g = (p & 0xff00) >> 8;
				int b = p & 0xff;
				data[i] = luminance(r, g, b);
			}
		}
		//For Black and White Image
		else if (type == BufferedImage.TYPE_BYTE_GRAY) {
			byte[] pixels = (byte[]) sourceImage.getData().getDataElements(0,
					0, width, height, null);
			for (int i = 0; i < picsize; i++) {
				data[i] = (pixels[i] & 0xff);
			}
		}
		//For Gray image with short datatype
		else if (type == BufferedImage.TYPE_USHORT_GRAY) {
			short[] pixels = (short[]) sourceImage.getData().getDataElements(0,
					0, width, height, null);
			for (int i = 0; i < picsize; i++) {
				data[i] = (pixels[i] & 0xffff) / 256;
			}
		}
		//For a Gray image with three each pixel
		else if (type == BufferedImage.TYPE_3BYTE_BGR) {
			byte[] pixels = (byte[]) sourceImage.getData().getDataElements(0,
					0, width, height, null);
			int offset = 0;
			for (int i = 0; i < picsize; i++) {
				int b = pixels[offset++] & 0xff;
				int g = pixels[offset++] & 0xff;
				int r = pixels[offset++] & 0xff;
				data[i] = luminance(r, g, b);
			}
		} else {
			throw new IllegalArgumentException("Unsupported image type: "
					+ type);
		}
	}

	/**
	 * Here we normalize the contrast of the image by making a histogram
	 */
	private static void normalizeContrast() {
		int[] histogram = new int[256];
		for (int i = 0; i < data.length; i++) {
			histogram[data[i]]++;
		}
		
		// All Histogram values are remap in 255 scale suxh that the distribution of all 
		// pixels become normalize in the image
		int[] remap = new int[256];
		int sum = 0;
		int j = 0;
		for (int i = 0; i < histogram.length; i++) {
			sum += histogram[i];
			int target = sum * 255 / picsize;
			for (int k = j + 1; k <= target; k++) {
				remap[k] = i;
			}
			j = target;
		}

		//Remapping is performed here
		for (int i = 0; i < data.length; i++) {
			data[i] = remap[data[i]];
		}
	}

	private static void writeEdges(int pixels[]) {

		if (edgesImage == null) {
			edgesImage = new BufferedImage(width, height,
					BufferedImage.TYPE_INT_ARGB);
		}
		edgesImage.getWritableTile(0, 0).setDataElements(0, 0, width, height,
				pixels);

	}

	public static void main(String args[]) throws Exception {

		/*
		 * 
		 * 
		 * Initializing Comm, world and size values here
		 */
		Comm.init(args);
		world = Comm.world();
		size = world.size();
		rank = world.rank();
		start_time= System.currentTimeMillis();
		
		// ...................... To read the image ...........................
		image = ImageIO.read(new File(args[0]));

		setLowThreshold(1f);
		setHighThreshold(2f);
		// apply it to an image

		ranges = new Range(0, height - 1).subranges(size);
		myrange = ranges[rank];
		mylb = myrange.lb();
		myub = myrange.ub();
		height = myub - mylb;

		Image img = image;
		int x = 0, y = 0;
		height = image.getHeight(null);
		width = image.getWidth(null);

		pixels = new int[width * height];

 
		PixelGrabber pg = new PixelGrabber(img, x, y, width, height, pixels, 0,
				width);
		try {
			pg.grabPixels();
		} catch (InterruptedException e) {
			System.err.println("interrupted waiting for pixels!");
			return;
		}
		if ((pg.getStatus() & ImageObserver.ABORT) != 0) {
			System.err.println("image fetch aborted or errored");
			return;
		}

	 
		// apply it to an image
		setSourceImage((BufferedImage) image);
		// detector.process();
/*
 * 
 * Process method is here;
 */
	 
		 
		picsize = width * height;
		//initialize array....
		if (data == null || picsize != data.length) {
			data = new int[picsize];
			magnitude = new int[picsize];

			xConv = new float[picsize];
			yConv = new float[picsize];
			xGradient = new float[picsize];
			yGradient = new float[picsize];
		}
		
		 
		readLuminance();
		
		if (contrastNormalized)
			normalizeContrast();
		
		
		
		computeGradients();
		
		
		
		int low = Math.round(lowThreshold * MAGNITUDE_SCALE);
		int high = Math.round(highThreshold * MAGNITUDE_SCALE);
 
		
		Arrays.fill(data, 0);

		int offset = 0;
		for ( y = 0; y < height; y++) {
			for ( x = 0; x < width; x++) {
				if (data[offset] == 0 && magnitude[offset] >= high) {
					follow(x, y, offset, low);
				}
				offset++;
			}
		}
		
		 
		
		for (int i = 0; i < picsize; i++) {
			data[i] = data[i] > 0 ? -1 : 0xff000000;
		}

	 
		IntegerBuf buf = IntegerBuf.buffer(data);

		world.reduce(0, buf, IntegerOp.SUM);
		
		
	/*
	 * Rank 0 performs the writing of the final output	
	 */
		if (rank == 0) {
			writeEdges(data);
			BufferedImage edges = getEdgesImage();
			end_time= System.currentTimeMillis();
			
			// Output written to the file
			File outputfile = new File("CannyoutputClu.png");
			ImageIO.write(edges, "png", outputfile);
			System.out.println("Time Taken: " + (end_time - start_time)+" msec");
		}

	}

}
