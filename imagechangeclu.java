/**
 * The Cluster Image Processing program
 * This program performs the image transformations like negative, 
 * sepia, flip, black and white image
 * 
 * The input image is sliced among the processors and then reduced after
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

import edu.rit.mp.IntegerBuf;
import edu.rit.pj.Comm;
import edu.rit.pj.reduction.IntegerOp;
import edu.rit.util.Range;

public class imagechangeclu
{
	static long time1,time2; 
	
	//Ranges for cluster processors
    static Comm world;
	static int size;
	static int rank;
	static Range[] ranges;
    static Range myrange;
    static int mylb;
    static int myub;

 	public static void main(String args[]) throws Exception
 	{
 		 Comm.init(args);
 		 world = Comm.world(); 
 		 size = world.size();
 		 rank = world.rank();
 		 	
		 Image image = null;
		 int height = 0, width = 0;
		 try 
		 {
		     // Read from a file
		     File file = new File(args[0]);
		     image = ImageIO.read(file);
		
		     height = image.getHeight(null);
		     width = image.getWidth(null);
		
		 } 
		 catch (IOException e) 
		 {
		 }
		
		 //Upper and lower bound calculations for processors
		 ranges = new Range (0, height-1).subranges (size);
		 myrange = ranges[rank];
		 mylb = myrange.lb();
		 myub = myrange.ub();
		 
		 Image img=image;
		 int x=0,y=0;
		 final int w=width;
		 final int h=height;
		 
		 final int[] pixels = new int[w * h];
			
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
		
		 time1=System.currentTimeMillis();
		
		final int[] negpixels = new int[w * h];
		final int [] blackpixels = new int[w * h];
		final int [] flippixels = new int[w * h];
		final int [] sepiapixels = new int[w * h];
		
		// loop local variable declaration
		for (int j = mylb; j <= myub; j++)
		{ 
			for (int i = 0; i < w; i++) 
		    {
				flippixels[j * w + i] = pixels[j * w + (w - i - 1)];
			 	int pixel =  pixels[j * w + i];
			 	int red   = (pixel >> 16) & 0xff;
				int green = (pixel >>  8) & 0xff;
				int blue  = (pixel      ) & 0xff;
				
				int blackwhite = (int) (red*0.3 + green*0.59+ blue*0.11);//red + green + blue;		
				blackpixels[j * w + i] = (blackwhite << 16) | (blackwhite <<  8) | (blackwhite);
				negpixels[j * w + i] = ((255-red) << 16) | ((255-green) <<  8) | (255-blue);
				
				int sepia = 20;
				int gray = (red + green + blue) / 3;
				red = green = blue = gray;
				red = red + (sepia * 2);
				green = green + sepia;

				if (red > 255) red=255;
				if (green > 255) green=255;
				if (blue > 255) blue=255;

				int sepiaIntensity = 20;
				// Darken blue color to increase sepia effect
				blue-= sepiaIntensity;

				// normalize if out of bounds
				if (blue < 0) blue=0;
				if (blue > 255) blue=255;
				
				sepiapixels[j * w + i] = ((red) << 16) | ((green) <<  8) | (blue);
		    }
		}
		
		//Buffer for the image to send
        IntegerBuf buf = IntegerBuf.buffer(negpixels); 
        IntegerBuf buf1 = IntegerBuf.buffer(flippixels); 
        IntegerBuf buf2 = IntegerBuf.buffer(blackpixels); 
        IntegerBuf buf3 = IntegerBuf.buffer(sepiapixels);
        
        //Reduce the results to the final matrix for image
        world.reduce(0,buf,IntegerOp.SUM);
        world.reduce(0,buf1,IntegerOp.SUM);
        world.reduce(0,buf2,IntegerOp.SUM);
        world.reduce(0,buf3,IntegerOp.SUM);
        
        if(rank == 0)
        {
        	BufferedImage pixelImage = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);    
    		pixelImage.setRGB(0, 0, w, h, blackpixels, 0, w);
    		    
    		File outputfile = new File("blackwhite.png");
    		ImageIO.write(pixelImage, "png", outputfile);
    		 
    		BufferedImage negpixelImage = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);    
    		negpixelImage.setRGB(0, 0, w, h, negpixels, 0, w);
    		
    	    // To write the Canny Edges to the image buffer...
    		File outputfile1 = new File("negative.png");
    		ImageIO.write(negpixelImage, "png", outputfile1);
    	    
    	    BufferedImage flippixelImage = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);    
    	    flippixelImage.setRGB(0, 0, w, h, flippixels, 0, w);

    		// retrieve image
    		File outputfile11 = new File("flip.png");
    		ImageIO.write(flippixelImage, "png", outputfile11);
    	    
    	    BufferedImage sepiapixelImage = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);    
    	    sepiapixelImage.setRGB(0, 0, w, h, sepiapixels, 0, w);

    		// retrieve image
    		File outputfile111 = new File("sepia.png");
    		ImageIO.write(sepiapixelImage, "png", outputfile111);
	    
	    	time2=System.currentTimeMillis();		
		 
			System.out.println((time2-time1)+" msec");
        }
	}
}
