package de.embl.cba.morphometry.algorithm;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import inra.ijpb.morphology.Strel;
import inra.ijpb.morphology.directional.DirectionalFilter.Operation;
import inra.ijpb.morphology.directional.OrientedLineStrelFactory;

public class DirectionalLineFilter
{

	private ImagePlus result;

	public void run( ImagePlus imagePlus, Operation operation, int lineLength, int nDirections )
	{
		ImageProcessor image = imagePlus.getProcessor();

		ImageStack resultStack = allocateResultImageStack( nDirections, image );

		OrientedLineStrelFactory strelFactory = new OrientedLineStrelFactory( lineLength );
		
		// Iterate over the set of directions
		for ( int i = 0; i < nDirections; i++)
		{
//			fireProgressChanged(this, i, nDirections);
			IJ.showProgress(i, nDirections );
			
			// Create the structuring element for current orientation
			double theta = ((double) i) * 180.0 / nDirections;
			Strel strel = strelFactory.createStrel(theta);

			// Apply oriented filter
			ImageProcessor filtered = operation.apply(image, strel);

			resultStack.setProcessor(filtered, i+1);
		}

		IJ.showProgress(1, 1);
		
		String name = imagePlus.getShortTitle() + "-orient";
		result = new ImagePlus(name, resultStack );
	}

	public ImageStack allocateResultImageStack( int nDirections, ImageProcessor image )
	{
		int sizeX = image.getWidth();
		int sizeY = image.getHeight();
		return ImageStack.create(sizeX, sizeY, nDirections, image.getBitDepth());
	}

	public ImagePlus getResult()
	{
		return result;
	}
}
