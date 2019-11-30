package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.Logger;
import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.plugin.*;
import ij.process.ColorProcessor;
import ij.process.LUT;
import loci.common.DebugTools;
import loci.formats.FormatException;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.widget.Button;
import loci.plugins.BF;

import java.awt.*;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Random;

@Plugin(type = Command.class, menuPath = "Plugins>EMBL>BD Processing" )
public class BDImageProcessingCommand implements Command
{
	@Parameter
	public LogService logService;

	@Parameter ( style = "directory", label = "Input Directory" )
	public File inputDirectory;

	@Parameter ( style = "directory", label = "Output Directory" )
	public File outputImageDirectory;

	@Parameter ( label = "Minimum File Size [kb]")
	public double minimumFileSizeKiloBytes = 10;

	@Parameter ( label = "Minimum Brightfield Intensity" )
	public double minBF = 0.08;

	@Parameter ( label = "Maximum Brightfield Intensity" )
	public double maxBF = 0.5;

	@Parameter ( label = "Minimum GFP Intensity" )
	public double minGFP = 0.08;

	@Parameter ( label = "Maximum GFP Intensity" )
	public double maxGFP = 1.0;

	@Parameter ( label = "Simple Overlay" )
	public boolean isSimpleOverlay = false;

	@Parameter ( label = "Preview Random Image", callback = "showRandomImage" )
	private Button showRandomImage;

	private String[] fileNames;
	private int numFiles;
	private ImagePlus inputImp;
	private ImagePlus outputImp;

	private ImagePlus brightFieldImp;
	private ImagePlus sideScatterImp;
	private ImagePlus forwardScatterImp;
	private ImagePlus gfpImp;
	private ImagePlus brightFieldGfpImp;

	public void run()
	{
		final long startMillis = System.currentTimeMillis();

		DebugTools.setRootLevel("OFF"); // Bio-Formats

		fetchFiles();

		createAndSaveImages();

		logFinished( startMillis );
	}

	private void createAndSaveImages()
	{
		for ( String fileName : fileNames )
		{
			if ( ! createOutputImage( fileName ) ) continue;
			saveImage( fileName, outputImp );
		}
	}

	private boolean checkFileSize( String filePath )
	{
		final File file = new File( filePath );
		final double fileSizeKiloBytes = getFileSizeKiloBytes( file );

		if ( fileSizeKiloBytes < minimumFileSizeKiloBytes )
		{
			IJ.log( "Skipped too small file: " + file.getName() );
			return false;
		}
		else
		{
			return true;
		}
	}

	private void showRandomImage()
	{
		if ( outputImp != null ) outputImp.close();

		fetchFiles();

		final String filePath = getRandomFilePath();

		if ( ! createOutputImage( filePath ) ) return;

		outputImp.show();
	}

	private String getRandomFilePath()
	{
		final Random random = new Random();
		return inputDirectory + File.separator + fileNames[ random.nextInt( numFiles ) ];
	}

	private void fetchFiles()
	{
		if ( fileNames != null ) return;

		final long startMillis = System.currentTimeMillis();
		fileNames = getValidFileNames();
		numFiles = fileNames.length;
		IJ.log( "Fetched file list in " + ( System.currentTimeMillis() - startMillis) + " ms; number of files: " + numFiles );
	}

	private static double getFileSizeKiloBytes( File file )
	{
		return file.length() / 1024.0 ;
	}

	private String[] getValidFileNames()
	{
		return inputDirectory.list( new FilenameFilter()
			{
				@Override
				public boolean accept( File dir, String name )
				{
					if (  ! name.contains( ".tif" ) ) return false;

					return true;
				}
			} );
	}

	private boolean createOutputImage( String filePath )
	{
		if ( ! checkFileSize( filePath ) ) return false;

		inputImp = tryOpenImage( filePath );

		inputImp = processImage( inputImp );

		extractChannels( inputImp );

		createRGBImages();

		outputImp = createOutputImp();

		outputImp.setTitle( new File( filePath ).getName() );

		return true;
	}

	private ImagePlus tryOpenImage( String filePath )
	{
		ImagePlus inputImp = null;
		try
		{
			inputImp = openImage( filePath );
		} catch ( FormatException e )
		{
			e.printStackTrace();
		} catch ( IOException e )
		{
			e.printStackTrace();
		}
		return inputImp;
	}

	private void extractChannels( ImagePlus inputImp )
	{
		brightFieldImp = getByteImagePlus( "bf", minBF, maxBF, 1, inputImp );
		sideScatterImp = getByteImagePlus( "ss", 0.08, 1.0, 2, inputImp );
		forwardScatterImp = getByteImagePlus( "fs", 0.08, 1.0, 4, inputImp );
		gfpImp = getByteImagePlus( "gfp", minGFP, maxGFP, 5, inputImp );
	}

	private ImagePlus openImage( String filePath ) throws FormatException, IOException
	{
		final ImagePlus[] imps = BF.openImagePlus( filePath );
		return imps[ 0 ];
	}

	private ImagePlus processImage( ImagePlus inputImp )
	{
		IJ.run(inputImp, "Scale...", "x=1.0 y=0.8 z=1.0 interpolation=Bilinear average process title=stack");
		IJ.run(inputImp, "Convolve...", "text1=[0 1.6 4 1.6 0\n] normalize stack");

		return inputImp;
	}

	private void showImages( String filePath, ImagePlus inputImp, ImagePlus outputImp )
	{
		final String fileName = new File( filePath ).getName();
		inputImp.setTitle( fileName );
		inputImp.show();
		outputImp.setTitle( fileName + "-output" );
		outputImp.show();
	}

	private void saveImage( String fileName, ImagePlus outputImp )
	{
		final String outputFileName = fileName.replace( ".tiff", "" );
		final String outputFilePath = outputImageDirectory + File.separator + outputFileName + ".jpg";
		new File( outputImageDirectory.toString() ).mkdirs();
		new FileSaver( outputImp ).saveAsJpeg( outputFilePath  );
	}

	private void createRGBImages()
	{
		brightFieldGfpImp = createBrightfieldGfpRGB( brightFieldImp, gfpImp );
		brightFieldImp = toRGB( brightFieldImp );
		sideScatterImp = toRGB( sideScatterImp );
		forwardScatterImp = toRGB( forwardScatterImp );
		gfpImp = toRGB( gfpImp );
	}

	private ImagePlus toRGB( ImagePlus imp )
	{
		return new ImagePlus( "", imp.getProcessor().convertToColorProcessor() );
	}

	private ImagePlus createOutputImp()
	{
		if ( isSimpleOverlay )
		{
			return brightFieldGfpImp;
		}
		else
		{
			// make montage
			final int width = brightFieldImp.getWidth();
			final int height = brightFieldImp.getHeight();

			int montageWidth = 3 * width;
			int montageHeight = 2 * height;

			final ImagePlus montageImp = new ImagePlus( "montage", new ColorProcessor( montageWidth, montageHeight ) );

			final StackInserter inserter = new StackInserter();
			inserter.insert( brightFieldImp, montageImp, 0, 0 );
			inserter.insert( gfpImp, montageImp, width, 0 );
			inserter.insert( brightFieldGfpImp, montageImp, 2 * width, 0 );
			inserter.insert( sideScatterImp, montageImp, 0, height );
			inserter.insert( forwardScatterImp, montageImp, width, height );

			return montageImp;
		}
	}

	private CompositeImage createBrightfieldGfpRGB( ImagePlus bfImp, ImagePlus gfpImp )
	{
		final CompositeImage bfGfpImp = ( CompositeImage ) RGBStackMerge.mergeChannels( new ImagePlus[]{ bfImp, gfpImp }, true );
		bfGfpImp.setC( 1 );
		bfGfpImp.setChannelLut( LUT.createLutFromColor( Color.WHITE ) );
		bfGfpImp.setDisplayRange( 0, 255 );
		bfGfpImp.setC( 2 );
		bfGfpImp.setChannelLut( LUT.createLutFromColor( Color.GREEN ) );
		bfGfpImp.setDisplayRange( 0, 255 );

		RGBStackConverter.convertToRGB( bfGfpImp );

		return bfGfpImp;
	}

	private static ImagePlus getByteImagePlus( String title, double minBF, double maxBF, int slice, ImagePlus imp )
	{
		Duplicator duplicator = new Duplicator();
		final ImagePlus impBf = duplicator.run( imp, slice, slice );
		impBf.getProcessor().setMinAndMax( minBF, maxBF );
		return new ImagePlus( title, impBf.getProcessor().convertToByteProcessor() );
	}

	private void logFinished( double startMillis )
	{
		IJ.log( "Processed " + numFiles + " in " + ( System.currentTimeMillis() - startMillis ) + " ms." );
		Logger.log( "Done!" );
	}


}
