package de.embl.cba.morphometry;

import ij.ImagePlus;
import ij.io.FileSaver;
import ij.measure.Calibration;
import loci.formats.FormatException;
import loci.plugins.in.ImagePlusReader;
import loci.plugins.in.ImportProcess;
import loci.plugins.in.ImporterOptions;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.view.Views;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;

public class ImageIO
{

	public static ImagePlus openWithBioFormats( String path )
	{
		try
		{
			ImporterOptions opts = new ImporterOptions();
			opts.setId( path );
			opts.setVirtual( true );

			ImportProcess process = new ImportProcess( opts );
			process.execute();

			ImagePlusReader impReader = new ImagePlusReader( process );

			ImagePlus[] imps = impReader.openImagePlus();
			return imps[ 0 ];
		}
		catch ( Exception e )
		{
			e.printStackTrace();
			return null;
		}

	}

	public static < T extends RealType< T > & NativeType< T > >
	void saveLabels(
			ArrayList< RandomAccessibleInterval< T > > labelings,
			Calibration calibration,
			String outputLabelingsPath )
	{
		final ImagePlus labelsImp = Utils.labelingsAsImagePlus( labelings );
		labelsImp.setCalibration( calibration );

		new FileSaver( labelsImp ).saveAsTiff( outputLabelingsPath );

		Logger.log( "Label images saved: " + outputLabelingsPath );
	}

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
