package de.embl.cba.morphometry.drosophila.registration;

import de.embl.cba.morphometry.Utils;
import ij.ImagePlus;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

public abstract class DrosophilaRegistrationUtils
{
	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > getChannelImages( ImagePlus imagePlus )
	{
		RandomAccessibleInterval< T > images = ImageJFunctions.wrap( imagePlus );

		int numChannels = imagePlus.getNChannels();

		if ( numChannels == 1 )
		{
			images = Views.addDimension( images, 0 ,0 );
		}
		else
		{
			images = Views.permute( images, Utils.imagePlusChannelDimension, 3 );
		}

		return images;
	}

	public static < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< T > getChannelImage( RandomAccessibleInterval< T > images, int channel )
	{
		RandomAccessibleInterval< T > rai = Views.hyperSlice( images, 3, channel );
		return rai;
	}
}
