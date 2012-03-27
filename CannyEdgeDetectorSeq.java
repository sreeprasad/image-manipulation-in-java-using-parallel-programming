/**
 * This program performs the Canny Edge Detection in a Sequential fashion.
 * It mainly detects the prominent edges in the image and display then in the output.
 * 
 * @author Divin Visariya
 * @author Sreeprasad Govindankutty
 * @author Varun Goyal
 * 
 */

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Arrays;

import javax.imageio.ImageIO;

/**
 * This class consits of methods which extracts the edges from the input image
 * 
 * The program takes input image as argument and outputs the image as 
 * "CannyoutputSeq.png"
 * 
 */

public class CannyEdgeDetectorSeq
{
	// To store the start and end time
	static long startTime,endTime;
	static int initX, maxX, initY, maxY;

	// Constants
	final static float GAUSSIAN_CUT_OFF = 0.005f;
	final static float MAGNITUDE_SCALE = 100F;
	final static float MAGNITUDE_LIMIT = 1000F;
	final static int MAGNITUDE_MAX = (int) (MAGNITUDE_SCALE * MAGNITUDE_LIMIT);

	// To store image properties
	static int height;
	static int width;
	static int picsize;
	static int[] pixelData;
	static int[] gradManitude;

	// Image buffers
	static BufferedImage sourceImage;
	static BufferedImage destEdgesImage;

	// To store the threshold values
	static float lowThreshold, highThreshold;
	static int low, high;

	// To store the gaussian kernel values
	static float gaussianKernelRadius;
	static int gaussianKernelWidth;

	// To store the gradient values of the image...
	static float[] xConv;
	static float[] yConv;
	static float[] xGradient;
	static float[] yGradient;

	/**
	 * The main program.
	 * 
	 * @param	args[]	Command line argument (ignored)
	 * 
	 */

	public static void main(String args[]) throws Exception
	{
		// ...................... To read the image ...........................
		sourceImage = ImageIO.read(new File(args[0]));
		
		// ............... To apply Canny Edge to the image ...................

		// To get picSize
		width = sourceImage.getWidth();
		height = sourceImage.getHeight();
		picsize = width * height;

		// To set the threshold values
		lowThreshold = 1f;
		highThreshold = 2f;

		// To set the Gaussian Kernel properties
		gaussianKernelRadius = 2f;
		gaussianKernelWidth = 16;

		// To initialize the arrays
		if (pixelData == null || picsize != pixelData.length) 
		{
			pixelData = new int[picsize];
			gradManitude = new int[picsize];

			xConv = new float[picsize];
			yConv = new float[picsize];

			xGradient = new float[picsize];
			yGradient = new float[picsize];
		}

		// To read the Luminance of the Source Image
		readLuminance();

		// To normalize the contrast of the image
		normalizeContrast();

		// To compute the start time
		startTime = System.currentTimeMillis();

		// To compute the gradients of the image
		computeGradients();

		// To compute the end time
		endTime = System.currentTimeMillis();

		// To perform Hysteresis
		Arrays.fill(pixelData, 0);
		low = Math.round(lowThreshold * MAGNITUDE_SCALE);
		high = Math.round( highThreshold * MAGNITUDE_SCALE);

		for (int y = 0, offset = 0; y < height; y++)
			for (int x = 0; x < width; x++, offset++) 
				if (pixelData[offset] == 0 && gradManitude[offset] >= high) 
					follow(x, y, offset, low);

		// To get Threshold edges
		for (int i = 0; i < picsize; i++) 
			pixelData[i] = (pixelData[i] > 0 ? -1 : 0xff000000);

		// To write the Canny Edges to the image buffer
		destEdgesImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		destEdgesImage.getWritableTile(0, 0).setDataElements(0, 0, width, height, pixelData);
		
		// retrieve image
		File outputfile = new File("CannyoutputSeq.png");
		ImageIO.write(destEdgesImage, "png", outputfile);

		// To print time taken
		System.out.println("Time Taken: " + (endTime - startTime)+" msec");
	}

	/**
	 * To compute the gradient values of the image.
	 * 
	 * @throws Exception
	 */

	static void computeGradients() throws Exception 
	{
		// To generate the Gaussian Convolution masks
		float kernel[] = new float[gaussianKernelWidth];
		float diffKernel[] = new float[gaussianKernelWidth];
		int kwidth;

		// Here the Gaussian mask is made using the formula
		for (kwidth = 0; kwidth < gaussianKernelWidth; kwidth++) 
		{
			float x = kwidth;
			float g1 = (float) Math.exp(-(x * x) / (2f * gaussianKernelRadius * gaussianKernelRadius));

			if (g1 <= GAUSSIAN_CUT_OFF && kwidth >= 2) 
				break;

			x = kwidth - 0.5f;
			float g2 = (float) Math.exp(-(x * x) / (2f * gaussianKernelRadius * gaussianKernelRadius));

			x = kwidth + 0.5f;
			float g3 = (float) Math.exp(-(x * x) / (2f * gaussianKernelRadius * gaussianKernelRadius));

			kernel[kwidth] = (g1 + g2 + g3) / 3f / (2f * (float) Math.PI * gaussianKernelRadius * gaussianKernelRadius);
			diffKernel[kwidth] = g3 - g2;
		}


		initX = kwidth - 1;
		maxX = width - (kwidth - 1);

		initY = width * (kwidth - 1);
		maxY = width * (height - (kwidth - 1));

		// To perform Convolution in x and y directions
		for (int x = initX; x < maxX; x++) 
		{
			for (int y = initY; y < maxY; y += width) 
			{
				int index = x + y;
				float sumX = pixelData[index] * kernel[0];
				float sumY = sumX;

				int yOffset = width;
				for(int xOffset=1; xOffset < kwidth; xOffset++) 
				{
					sumY += kernel[xOffset] * (pixelData[index - yOffset] + pixelData[index + yOffset]);
					sumX += kernel[xOffset] * (pixelData[index - xOffset] + pixelData[index + xOffset]);
					yOffset += width;	
				}

				xConv[index] = sumX;
				yConv[index] = sumY;
			}

		}

		// Here we calculate the xGradient values
		for (int x = initX; x < maxX; x++) 
		{
			for (int y = initY; y < maxY; y += width) 
			{
				float sum = 0f;
				int index = x + y;
				for (int i = 1; i < kwidth; i++)
					sum += diffKernel[i] * (yConv[index - i] - yConv[index + i]);

				xGradient[index] = sum;
			}
		}

		// Here we calculate the yGradient values
		for (int x = kwidth; x < width - kwidth; x++) 
		{
			for (int y = initY; y < maxY; y += width) 
			{
				float sum = 0.0f;
				int index = x + y;
				int yOffset = width;
				for (int i = 1; i < kwidth; i++) {
					sum += diffKernel[i] * (xConv[index - yOffset] - xConv[index + yOffset]);
					yOffset += width;
				}
				yGradient[index] = sum;
			}
		}


		initX = kwidth;
		maxX = width - kwidth;

		initY = width * kwidth;
		maxY = width * (height - kwidth);


		for (int x = initX; x < maxX; x++) 
		{
			for (int y = initY; y < maxY; y += width) 
			{
				//variables with all 8 directions stored around a pixel index
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

				float gradMag = (float) Math.sqrt((xGrad*xGrad) + (yGrad*yGrad));
				
				// Here we are finding the magnitude of the pixels in north, east, west, south and
				// north-east, north-west, south-east and sounth-west directions.
				float nMag = (float) Math.sqrt(
						(xGradient[indexN]*xGradient[indexN])+
						(yGradient[indexN]*yGradient[indexN])); 

				float sMag =(float) Math.sqrt(
						(xGradient[indexS]*xGradient[indexS])+
						(yGradient[indexS]*yGradient[indexS]));

				float wMag =(float) Math.sqrt(
						(xGradient[indexW]*xGradient[indexW])+
						(yGradient[indexW]*yGradient[indexW])); 

				float eMag =(float) Math.sqrt(
						(xGradient[indexE]*xGradient[indexE])+
						(yGradient[indexE]*yGradient[indexE])); 

				float neMag = (float) Math.sqrt(
						(xGradient[indexNE]*xGradient[indexNE])+
						(yGradient[indexNE]*yGradient[indexNE])); 

				float seMag = (float) Math.sqrt(
						(xGradient[indexSE]*xGradient[indexSE])+
						(yGradient[indexSE]*yGradient[indexSE])); 

				float swMag =(float) Math.sqrt(
						(xGradient[indexSW]*xGradient[indexSW])+
						(yGradient[indexSW]*yGradient[indexSW])); 

				float nwMag = (float) Math.sqrt(
						(xGradient[indexNW]*xGradient[indexNW])+
						(yGradient[indexNW]*yGradient[indexNW])); 




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

				if (xGrad * yGrad <= (float) 0 /*(1)*/
						? Math.abs(xGrad) >= Math.abs(yGrad) /*(2)*/
								? (tmp = Math.abs(xGrad * gradMag)) >= Math.abs(yGrad * neMag - (xGrad + yGrad) * eMag) /*(3)*/
										&& tmp > Math.abs(yGrad * swMag - (xGrad + yGrad) * wMag) /*(4)*/
										: (tmp = Math.abs(yGrad * gradMag)) >= Math.abs(xGrad * neMag - (yGrad + xGrad) * nMag) /*(3)*/
										&& tmp > Math.abs(xGrad * swMag - (yGrad + xGrad) * sMag) /*(4)*/
										: Math.abs(xGrad) >= Math.abs(yGrad) /*(2)*/
										? (tmp = Math.abs(xGrad * gradMag)) >= Math.abs(yGrad * seMag + (xGrad - yGrad) * eMag) /*(3)*/
												&& tmp > Math.abs(yGrad * nwMag + (xGrad - yGrad) * wMag) /*(4)*/
												: (tmp = Math.abs(yGrad * gradMag)) >= Math.abs(xGrad * seMag + (yGrad - xGrad) * sMag) /*(3)*/
												&& tmp > Math.abs(xGrad * nwMag + (yGrad - xGrad) * nMag) /*(4)*/
				)

				{
					gradManitude[index] = gradMag >= MAGNITUDE_LIMIT ? MAGNITUDE_MAX : (int) (MAGNITUDE_SCALE * gradMag);
				} 
				else 
				{
					gradManitude[index] = 0;
				}
			}
		}
	}

	// Here we mainly follow the lines according to the threashold values
	// All pixels with magnitude above the treshhold are considered
	static void follow(int x1, int y1, int i1, int threshold) 
	{
		int x0 = x1 == 0 ? x1 : x1 - 1;
		int x2 = x1 == width - 1 ? x1 : x1 + 1;
		int y0 = y1 == 0 ? y1 : y1 - 1;
		int y2 = y1 == height -1 ? y1 : y1 + 1;

		pixelData[i1] = gradManitude[i1];
		for (int x = x0; x <= x2; x++) 
		{
			for (int y = y0; y <= y2; y++) 
			{
				int i2 = x + y * width;
				if ((y != y1 || x != x1) && pixelData[i2] == 0 && gradManitude[i2] >= threshold) 
				{
					follow(x, y, i2, threshold);
					return;
				}
			}
		}
	}

	/**
	 * This methods calculated the luminance of every pixel of the image.
	 * 
	 * The image are selected according to their types and luminance is 
	 * performed accordingly
	 * 
	 */
	static void readLuminance() 
	{
		int type = sourceImage.getType();

		//For RGB image
		if (type == BufferedImage.TYPE_INT_RGB || type == BufferedImage.TYPE_INT_ARGB) 
		{
			int[] pixels = (int[]) sourceImage.getData().getDataElements(0, 0, width, height, null);		
			for (int i = 0; i < picsize; i++) 
			{
				int p = pixels[i];
				int r = ((p & 0xff0000) >> 16);
				int g = ((p & 0xff00) >> 8);
				int b = (p & 0xff);
				pixelData[i] = Math.round(0.299f * r + 0.587f * g + 0.114f * b);
			} 
		} 
		//For Black and White Image
		else if (type == BufferedImage.TYPE_BYTE_GRAY) 
		{
			byte[] pixels = (byte[]) sourceImage.getData().getDataElements(0, 0, width, height, null);
			for (int i = 0; i < picsize; i++) 
			{
				pixelData[i] = (pixels[i] & 0xff);
			}
		}
		//For Gray image with short datatype
		else if (type == BufferedImage.TYPE_USHORT_GRAY) 
		{
			short[] pixels = (short[]) sourceImage.getData().getDataElements(0, 0, width, height, null);
			for (int i = 0; i < picsize; i++) 
			{
				pixelData[i] = (pixels[i] & 0xffff) / 256;
			}
		}
		//For a Gray image with three each pixel
		else if (type == BufferedImage.TYPE_3BYTE_BGR) 
		{
			byte[] pixels = (byte[]) sourceImage.getData().getDataElements(0, 0, width, height, null);
			int offset = 0;
			for (int i = 0; i < picsize; i++) 
			{
				int b = pixels[offset++] & 0xff;
				int g = pixels[offset++] & 0xff;
				int r = pixels[offset++] & 0xff;
				pixelData[i] = Math.round(0.299f * r + 0.587f * g + 0.114f * b);
			}
		}
		else 
		{
			throw new IllegalArgumentException("Unsupported image type: " + type);
		}
	}

	/**
	 * Here we normalize the contrast of the image by making a histogram
	 */
	static void normalizeContrast() 
	{
		int[] histogram = new int[256];
		for (int i = 0; i < pixelData.length; i++) 
		{
			histogram[pixelData[i]]++;
		}

		// All Histogram values are remap in 255 scale suxh that the distribution of all 
		// pixels become normalize in the image
		int[] remap = new int[256];
		int sum = 0;
		int j = 0;
		for (int i = 0; i < histogram.length; i++) 
		{
			sum += histogram[i];
			int target = sum*255/picsize;
			for (int k = j+1; k <=target; k++) 
			{
				remap[k] = i;
			}
			j = target;
		}

		//Remapping is performed here
		for (int i = 0; i < pixelData.length; i++) 
			pixelData[i] = remap[pixelData[i]];

	}

}
