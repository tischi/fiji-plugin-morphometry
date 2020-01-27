package de.embl.cba.morphometry.fccf;

import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackConverter;
import ij.plugin.RGBStackMerge;
import ij.plugin.StackInserter;
import ij.process.ColorProcessor;
import ij.process.LUT;
import loci.formats.FormatException;
import loci.plugins.BF;

import java.awt.*;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public abstract class FCCF
{
	public static final String BRIGHTFIELD_AND_GFP = "BrightfieldAndGFP";
	public static final String GREEN_FLUORESCENCE = "GreenFluorescence";
	public static final String FOREWARD_SCATTER = "ForewardScatter";
	public static final String SIDE_SCATTER = "SideScatter";
	public static final String BRIGHTFIELD = "BrightField";

	public static final String[] imageNames = new String[]{
			BRIGHTFIELD, SIDE_SCATTER, FOREWARD_SCATTER, GREEN_FLUORESCENCE
	};
	public static final String VIEW_RAW = "Raw";
	public static final String VIEW_PROCESSED_MONTAGE = "Processed Montage";
	public static final String VIEW_PROCESSED_BF_GF_OVERLAY = "Processed BrightField GreenFluo Overlay";
	public static final String VIEW_RAW_AND_MONTAGE = "Raw and Processed Montage";

	public static HashMap< String, Integer > getNameToSlice()
	{
		final HashMap< String, Integer > nameToSlice = new HashMap<>();

		nameToSlice.put( BRIGHTFIELD, 1 );
		nameToSlice.put( GREEN_FLUORESCENCE, 5 );
		nameToSlice.put( FOREWARD_SCATTER, 4 );
		nameToSlice.put( SIDE_SCATTER, 2 );

		return nameToSlice;
	}

	public static ImagePlus createProcessedImage(
			String filePath,
			Map< String, double[] > nameToRange,
			Map< String, Integer > nameToSlice,
			String viewingModality )
	{
		ImagePlus inputImp = tryOpenImage( filePath );

		if ( viewingModality.equals( FCCF.VIEW_RAW ) ) return inputImp;

		inputImp = processImage( inputImp );

		final Map< String, ImagePlus > nameToImp = extractChannels( inputImp, nameToRange, nameToSlice );

		convertImagesToRGB( nameToImp );

		ImagePlus outputImp = createOutputImp( nameToImp, viewingModality );

		outputImp.setTitle( new File( filePath ).getName() );

		return outputImp;
	}

	public static ImagePlus tryOpenImage( String filePath )
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

	public static HashMap< String, ImagePlus > extractChannels(
			ImagePlus inputImp,
			Map< String, double[] > nameToRange,
			Map< String, Integer > nameToSlice )
	{
		final HashMap< String, ImagePlus > nameToImagePlus = new HashMap<>();

		for( String name : imageNames )
			putImagePlus( inputImp, nameToImagePlus, name, nameToRange, nameToSlice );

		return nameToImagePlus;
	}

	public static void putImagePlus(
			ImagePlus inputImp,
			HashMap< String, ImagePlus > nameToImp,
			String name,
			Map< String, double[] > nameToRange,
			Map< String, Integer > nameToSlice )
	{
		final ImagePlus byteImagePlus =
				getByteImagePlus(
						name,
						nameToRange.get( name )[ 0 ],
						nameToRange.get( name )[ 1 ],
						nameToSlice.get( name ),
						inputImp );

		nameToImp.put( name, byteImagePlus );
	}

	public static ImagePlus openImage( String filePath ) throws FormatException, IOException
	{
		final ImagePlus[] imps = BF.openImagePlus( filePath );
		return imps[ 0 ];
	}

	public static ImagePlus processImage( ImagePlus inputImp )
	{
		IJ.run(inputImp, "Scale...", "x=1.0 y=0.8 z=1.0 interpolation=Bilinear average process title=stack");
		IJ.run(inputImp, "Convolve...", "text1=[0 1.6 4 1.6 0\n] normalize stack");

		return inputImp;
	}

	public static void convertImagesToRGB( Map< String, ImagePlus > nameToImp )
	{
		final CompositeImage brightfieldGfpRGB = createBrightfieldGfpRGB( nameToImp.get( BRIGHTFIELD ), nameToImp.get( GREEN_FLUORESCENCE ) );
		nameToImp.put( BRIGHTFIELD_AND_GFP, brightfieldGfpRGB );

		for ( String name : imageNames )
			nameToImp.put( name, toRGB( nameToImp.get( name ) ) );

	}

	public static ImagePlus toRGB( ImagePlus imp )
	{
		return new ImagePlus( "", imp.getProcessor().convertToColorProcessor() );
	}

	public static ImagePlus createOutputImp( final Map< String, ImagePlus > nameToImp, String viewingModality )
	{
		if ( viewingModality.equals( FCCF.VIEW_PROCESSED_BF_GF_OVERLAY ) )
		{
			return nameToImp.get( BRIGHTFIELD_AND_GFP );
		}
		else if ( viewingModality.equals( FCCF.VIEW_PROCESSED_MONTAGE ) )
		{
			// make montage
			final int width = nameToImp.get( BRIGHTFIELD ).getWidth();
			final int height = nameToImp.get( BRIGHTFIELD ).getHeight();

			int montageWidth = 3 * width;
			int montageHeight = 2 * height;

			final ImagePlus montageImp = new ImagePlus( "montage", new ColorProcessor( montageWidth, montageHeight ) );

			final StackInserter inserter = new StackInserter();
			inserter.insert( nameToImp.get( BRIGHTFIELD ), montageImp, 0, 0 );
			inserter.insert( nameToImp.get( GREEN_FLUORESCENCE ), montageImp, width, 0 );
			inserter.insert( nameToImp.get( BRIGHTFIELD_AND_GFP ), montageImp, 2 * width, 0 );
			inserter.insert( nameToImp.get( SIDE_SCATTER ), montageImp, 0, height );
			inserter.insert( nameToImp.get( FOREWARD_SCATTER ), montageImp, width, height );

			return montageImp;
		}
		else
		{
			throw new UnsupportedOperationException( "Viewing modality not supported: " + viewingModality );
		}
	}

	public static CompositeImage createBrightfieldGfpRGB( ImagePlus bfImp, ImagePlus gfpImp )
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

	public class IntensityRanges
	{
		double minBF;
		double maxBF;
		double minGFP;
		double maxGFP;
	}

	public static boolean checkFile( String filePath, double minimumFileSizeKiloBytes, double maximumFileSizeKiloBytes )
	{
		final File file = new File( filePath );

		if ( ! file.exists() )
		{
			throw new UnsupportedOperationException( "File does not exist: " + file );
		}

		final double fileSizeKiloBytes = getFileSizeKiloBytes( file );

		if ( fileSizeKiloBytes < minimumFileSizeKiloBytes )
		{
			IJ.log( "Skipped too small file: " + file.getName() + "; size [kB]: " + fileSizeKiloBytes);
			return false;
		}
		else if ( fileSizeKiloBytes > maximumFileSizeKiloBytes )
		{
			IJ.log( "Skipped too large file: " + file.getName() + "; size [kB]: " + fileSizeKiloBytes);
			return false;
		}
		else
		{
			return true;
		}
	}

	public static double getFileSizeKiloBytes( File file )
	{
		return file.length() / 1024.0 ;
	}

	public static ImagePlus getByteImagePlus( String title, double minBF, double maxBF, int slice, ImagePlus imp )
	{
		Duplicator duplicator = new Duplicator();
		final ImagePlus impBf = duplicator.run( imp, slice, slice );
		impBf.getProcessor().setMinAndMax( minBF, maxBF );
		return new ImagePlus( title, impBf.getProcessor().convertToByteProcessor() );
	}

	public static String[] getValidFileNames( File inputDirectory )
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
}
