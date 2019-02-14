package de.embl.cba.morphometry;

import ij.ImagePlus;
import ij.io.FileSaver;
import loci.formats.FormatException;
import loci.plugins.in.ImagePlusReader;
import loci.plugins.in.ImportProcess;
import loci.plugins.in.ImporterOptions;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;

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

	public static void saveImages( String inputPath, ArrayList< ImagePlus > imps )
	{
		for ( ImagePlus imp : imps )
		{
			final String outputPath = inputPath + "-" + imp.getTitle() + ".tif";
			FileSaver fileSaver = new FileSaver( imp );
			fileSaver.saveAsTiff( outputPath );
		}
	}

	public static < T extends RealType< T > & NativeType< T > >
	void saveLabels( ArrayList< RandomAccessibleInterval< T > > labelings, String outputLabelingsPath )
	{
		final ImagePlus labelsImp = Utils.labelingsAsImagePlus( labelings );

		new FileSaver( labelsImp ).saveAsTiff( outputLabelingsPath );

		Logger.log( "Label images saved: " + outputLabelingsPath );
	}
}
