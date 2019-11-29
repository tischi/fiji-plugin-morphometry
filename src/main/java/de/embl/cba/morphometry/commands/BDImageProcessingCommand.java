package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.spindle.SpindleMorphometrySettings;
import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.plugin.*;
import ij.process.ColorProcessor;
import ij.process.LUT;
import loci.formats.FormatException;
import net.imagej.ops.OpService;
import net.imglib2.type.numeric.RealType;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import loci.plugins.BF;

import java.awt.*;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;

@Plugin(type = Command.class, menuPath = "Plugins>EMBL>BD Processing" )
public class BDImageProcessingCommand< R extends RealType< R > > implements Command
{
	public SpindleMorphometrySettings settings = new SpindleMorphometrySettings();

	@Parameter
	public OpService opService;

	@Parameter ( style = "directory", label = "Input Directory" )
	public File inputImageDirectory;

	@Parameter ( style = "directory", label = "Output Directory" )
	public File outputImageDirectory;

	@Parameter ( label = "Minimum BF" )
	public double minBF = 0.08;

	@Parameter ( label = "Maximum BF" )
	public double maxBF = 0.5;

	@Parameter ( label = "Minimum GFP" )
	public double minGFP = 0.08;

	@Parameter ( label = "Maximum GFP" )
	public double maxGFP = 1.0;

	@Parameter ( label = "Simple Overlay" )
	public boolean isSimpleOverlay;
	private ImagePlus bfImp;
	private ImagePlus ssImp;
	private ImagePlus fsImp;
	private ImagePlus gfpImp;
	private ImagePlus bfGfpImp;

	public void run()
	{
		final String[] fileNames = inputImageDirectory.list( new FilenameFilter()
		{
			@Override
			public boolean accept( File dir, String name )
			{
				if ( name.contains( ".tif" ) ) return true;
				else return false;
			}
		} );

		for ( String fileName : fileNames )
		{
			tryProcessFile( fileName );
			break;
		}
	}

	private void tryProcessFile( String fileName )
	{
		try
		{
			processFile( inputImageDirectory + File.separator + fileName );
		} catch ( IOException e )
		{
			e.printStackTrace();
		} catch ( FormatException e )
		{
			e.printStackTrace();
		}
	}

	private void processFile( String filePath ) throws IOException, FormatException
	{
		final ImagePlus[] imps = BF.openImagePlus( filePath );
		final ImagePlus imp = imps[ 0 ];

		IJ.run(imp, "Scale...", "x=1.0 y=0.8 z=1.0 interpolation=Bilinear average process title=stack");
		IJ.run(imp, "Convolve...", "text1=[0 1.6 4 1.6 0\n] normalize stack");

		bfImp = getByteImagePlus( "bf", minBF, maxBF, 1, imp );
		ssImp = getByteImagePlus( "ss", 0.08, 1.0, 2, imp );
		fsImp = getByteImagePlus( "fs", 0.08, 1.0, 4, imp );
		gfpImp = getByteImagePlus( "fitc", minGFP, maxGFP, 5, imp );

		bfGfpImp = getBfGfpRGB( bfImp, gfpImp );

		convertImagesToRGB();

		final ImagePlus outputImp = createOutputImp();

		save( filePath, outputImp );

		outputImp.show();

		logEnd();
	}

	private void save( String filePath, ImagePlus outputImp )
	{
		final String outputFileName = new File( filePath ).getName().replace( ".tiff", "" );
		final String outputFilePath = outputImageDirectory + File.separator + outputFileName + ".jpg";
		new FileSaver( outputImp ).saveAsJpeg( outputFilePath  );
	}

	private void convertImagesToRGB()
	{
		bfImp = toRGB( bfImp );
		ssImp = toRGB( ssImp );
		fsImp = toRGB( fsImp );
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
			return bfGfpImp;
		}
		else
		{
			// make montage
			final int width = bfImp.getWidth();
			final int height = bfImp.getHeight();

			int montageWidth = 3 * width;
			int montageHeight = 2 * height;

			final ImagePlus montageImp = new ImagePlus( "montage", new ColorProcessor( montageWidth, montageHeight ) );

			final StackInserter inserter = new StackInserter();
			inserter.insert( bfImp, montageImp, 0, 0 );
			inserter.insert( gfpImp, montageImp, width, 0 );
			inserter.insert( bfGfpImp, montageImp, 2 * width, 0 );
			inserter.insert( ssImp, montageImp, 0, height );
			inserter.insert( fsImp, montageImp, width, height );

			return montageImp;
		}
	}

	private CompositeImage getBfGfpRGB( ImagePlus bfImp, ImagePlus gfpImp )
	{
		final CompositeImage bfGfpImp = ( CompositeImage ) RGBStackMerge.mergeChannels( new ImagePlus[]{ bfImp, gfpImp }, true );
		bfGfpImp.setC( 1 );
		bfGfpImp.setChannelLut( LUT.createLutFromColor( Color.WHITE ) );
		bfGfpImp.setDisplayRange( 0, 255 );
		bfGfpImp.setC( 2 );
		bfGfpImp.setChannelLut( LUT.createLutFromColor( Color.GREEN ) );
		bfGfpImp.setDisplayRange( 0, 255 );

		bfGfpImp.show();

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

	private void logEnd()
	{
		Logger.log( "Done!" );
	}


}
