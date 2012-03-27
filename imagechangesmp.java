/**
 * The Parallel Image Processing program for SMP computers
 * This program performs the image transformations like negative, 
 * sepia, flip, black and white image
 * 
 * The input image is sliced among the smp cores and then reduced after
 * the manipulation is done
 * 
 * @author Divin Visariya
 * @author Sreeprasad Govindankutty
 * @author Varun Goyal
 */

import java.awt.Image;
import java.awt.image.*;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import edu.rit.pj.Comm;
import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

public class imagechangesmp
{
	static long time1,time2; 
	static BufferedImage destEdgesImage;
	

	/**
	 * The main method
	 * 
	 * @param args (input file is taken as argument)
	 * @throws IOException
	 */
 	public static void main(String args[]) throws Exception
 	{
 		//Scheduler for multiple processors
 		Comm.init(args);
	
		 Image image = null;
		 int height = 0, width = 0;
		 try 
		 {
		     // Read from a file
		     File file = new File(args[0]);
		     image = ImageIO.read(file);
		     
		     //Width and Height of the Image		
		     height = image.getHeight(null);
		     width = image.getWidth(null);
		 } 
		 catch (IOException e) 
		 {
		 }
		
		 Image img=image;
		 int x=0,y=0;
		 final int w=width;
		 final int h=height;
		
		//Array for storing pixels
		 final int[] pixels = new int[w * h];
			
		//Pixelgrabber grabs the pixels from the image to the array
		PixelGrabber pg = new PixelGrabber(img, x, y, w, h, pixels, 0, w);
		try 
		{
		    pg.grabPixels();
		} 
		catch (InterruptedException e) 
		{
		    System.err.println("interrupted waiting for pixels!");
		    return;
		}
		if ((pg.getStatus() & ImageObserver.ABORT) != 0) 
		{
		    System.err.println("image fetch aborted or errored");
		    return;
		}
		
		//Start timer		
		time1=System.currentTimeMillis();
		
		//Arrays initialized for different output
		final int [] negpixels = new int[w * h];
		final int [] blackpixels = new int[w * h];
		final int [] flippixels = new int[w * h];
		final int [] sepiapixels = new int[w * h];
		
		//Parallel Job for executing the image transformation in parallel
		new ParallelTeam().execute(new ParallelRegion() 
		{
			public void run() 
			{
				try 
				{
					execute(0, h - 1, new IntegerForLoop() 
					{
						// per thread variable declaration 	
						public void run(int first, int last) 
						{
							// loop local variable declaration
							for (int j = first; j <= last; j++)
							{ 
								for (int i = 0; i < w; i++) 
							    {
									//Copying the pixels accordingly to get the flip version
									 flippixels[j * w + i] = pixels[j * w + (w - i - 1)];
									 
									 //Retrieving the Red, Green and Blue section of every pixel 
									 int pixel =  pixels[j * w + i];
									 int red   = (pixel >> 16) & 0xff;
									 int green = (pixel >>  8) & 0xff;
									 int blue  = (pixel) & 0xff;
									
									 //The calculations on RGB values to get the black and white equivalent
									 int blackwhite = (int) (red *0.3 + green * 0.59+ blue * 0.11);		
									 blackpixels[j * w + i] = (blackwhite << 16) | (blackwhite <<  8) | (blackwhite);
									 
									 //Calculating the negative version of the pixel by inverting it
									 negpixels[j * w + i] = ((255-red) << 16) | ((255-green) <<  8) | (255-blue);
									
									 //The amount of sepia transformation to be done
									 int sepia = 20;
									 
									 //Calculation to get the sepia form from RGB values 
									 int gray = (red + green + blue) / 3;
									 red = green = blue = gray;
									 red = red + (sepia * 2);
									 green = green + sepia;

									 //Rounding of to 255
									 if (red > 255) red=255;
									 if (green > 255) green=255;
									 if (blue > 255) blue=255;

									 //Variable to increase the sepia effect
									 int sepiaIntensity = 20;
									 
									 // Darken blue color to increase sepia effect
									 blue-= sepiaIntensity;

									 // normalize if out of bounds
									 if (blue < 0) blue=0;
									 if (blue > 255) blue=255;
									 
									 sepiapixels[j * w + i] = ((red) << 16) | ((green) <<  8) | (blue);
							    }
							}

						}
					});
				} 
				catch (Exception e) 
				{
					e.printStackTrace();
				}
			}
		});
	
		//End Timer
		time2=System.currentTimeMillis();		
		 
		//Time Taken 
		System.out.println((time2-time1)+" msec");
		
		//Store all the ouput in the image files
		 
		 BufferedImage pixelImage = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);    
		 pixelImage.setRGB(0, 0, w, h, blackpixels, 0, w);
		 
		 //Black and white format output
		 File outputfile = new File("blackwhite.png");
		 ImageIO.write(pixelImage, "png", outputfile);
		 
		 BufferedImage negpixelImage = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);    
		 negpixelImage.setRGB(0, 0, w, h, negpixels, 0, w);
		
		 //Inverted(negative) format output
		 File outputfile1 = new File("negative.png");
		 ImageIO.write(negpixelImage, "png", outputfile1);
	    
		 BufferedImage flippixelImage = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);    
		 flippixelImage.setRGB(0, 0, w, h, flippixels, 0, w);

		 //Flip horizontally output
		 File outputfile11 = new File("flip.png");
		 ImageIO.write(flippixelImage, "png", outputfile11);
	    
		 BufferedImage sepiapixelImage = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);    
		 sepiapixelImage.setRGB(0, 0, w, h, sepiapixels, 0, w);

		 //Sepia format output
		 File outputfile111 = new File("sepia.png");
		 ImageIO.write(sepiapixelImage, "png", outputfile111);
	    
	}
}
