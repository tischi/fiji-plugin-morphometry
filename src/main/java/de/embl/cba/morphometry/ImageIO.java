package de.embl.cba.morphometry;

import ij.ImagePlus;
import ij.io.FileSaver;
import loci.formats.FormatException;
import loci.plugins.in.ImagePlusReader;
import loci.plugins.in.ImportProcess;
import loci.plugins.in.ImporterOptions;

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
}
