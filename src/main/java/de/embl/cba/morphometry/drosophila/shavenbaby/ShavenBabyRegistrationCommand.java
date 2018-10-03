package de.embl.cba.morphometry.drosophila.shavenbaby;

import bdv.util.*;
import de.embl.cba.morphometry.Projection;
import de.embl.cba.morphometry.refractiveindexmismatch.RefractiveIndexMismatchCorrectionSettings;
import de.embl.cba.morphometry.refractiveindexmismatch.RefractiveIndexMismatchCorrections;
import de.embl.cba.morphometry.Transforms;
import de.embl.cba.morphometry.Utils;
import ij.ImagePlus;
import ij.io.FileSaver;
import net.imagej.DatasetService;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.io.File;
import java.util.ArrayList;

import static de.embl.cba.morphometry.Constants.*;
import static de.embl.cba.morphometry.ImageIO.openWithBioFormats;


@Plugin(type = Command.class, menuPath = "Plugins>Registration>EMBL>Drosophila Shavenbaby" )
public class ShavenBabyRegistrationCommand <T extends RealType<T> & NativeType< T > > implements Command
{
	@Parameter
	public UIService uiService;

	@Parameter
	public DatasetService datasetService;

	@Parameter
	public LogService logService;

	@Parameter
	public OpService opService;

	@Parameter
	public StatusService statusService;

	@Parameter( required = false )
	public ImagePlus imagePlus;

	ShavenBabyRegistrationSettings settings = new ShavenBabyRegistrationSettings();

	public static final String FROM_DIRECTORY = "From directory";
	public static final String CURRENT_IMAGE = "Current image";

	@Parameter( choices = { FROM_DIRECTORY, CURRENT_IMAGE })
	public String inputModality = FROM_DIRECTORY;

	@Parameter( style = "directory" )
	public File inputDirectory;

	@Parameter
	public String fileNameEndsWith = ".czi,.lsm";

	@Parameter
	public boolean showIntermediateResults = settings.showIntermediateResults;

	@Parameter( choices = {
			ShavenBabyRegistrationSettings.CENTROID_SHAPE,
			ShavenBabyRegistrationSettings.PROJECTION_SHAPE,
			ShavenBabyRegistrationSettings.AMNIOSEROSA})
	public String rollAngleComputationMethod = settings.rollAngleComputationMethod;

	@Parameter
	public double registrationResolution = settings.registrationResolution;

	@Parameter
	public double outputResolution = settings.outputResolution;

//	@Parameter
//	public double thresholdInUnitsOfBackgroundPeakHalfWidth = settings.thresholdInUnitsOfBackgroundPeakHalfWidth;

//	@Parameter
//	public double refractiveIndexAxialCalibrationCorrectionFactor = settings.refractiveIndexAxialCalibrationCorrectionFactor;

	@Parameter
	public double refractiveIndexIntensityCorrectionDecayLength = settings.refractiveIndexIntensityCorrectionDecayLength;

//	@Parameter
//	public double watershedSeedsGlobalDistanceThreshold = settings.watershedSeedsGlobalDistanceThreshold;


	public void run()
	{
		setSettingsFromUI();

		final ShavenBabyRegistration registration = new ShavenBabyRegistration( settings, opService );


		if ( inputModality.equals( CURRENT_IMAGE ) && imagePlus != null )
		{
//			RandomAccessibleInterval< T > transformed = alignAndMaskImages( imagePlus, registration );
//			showWithBdv( transformed, "registered" );
//			ImageJFunctions.show( Views.permute( transformed, 2, 3 ) );
		}


		if ( inputModality.equals( FROM_DIRECTORY ) )
		{
			final File directory = inputDirectory;//uiService.chooseFile( null, FileWidget.DIRECTORY_STYLE );
			String[] files = directory.list();

			for( String file : files )
			{
				if ( acceptFile( fileNameEndsWith, file ) )
				{
					// Open
					final String inputPath = directory + "/" + file;
					Utils.log( " " );
					Utils.log( "Reading: " + inputPath + "..." );
					final ImagePlus inputImagePlus = openWithBioFormats( inputPath );

					if ( inputImagePlus == null )
					{
						logService.error( "Error opening file: " + inputPath );
						continue;
					}

					RandomAccessibleInterval< T > registeredImages = alignAndMaskImages( inputImagePlus, registration );

					// Save watershed
					RandomAccessibleInterval< T > watershed = (RandomAccessibleInterval) registration.getWatershedLabelImg();
					new FileSaver( ImageJFunctions.wrap( watershed, "" ) ).saveAsTiff( inputPath + "-watershed.tif" );

					if ( registeredImages == null )
					{
						Utils.log( "ERROR: Could not find central embryo" );
						continue;
					}

					if ( settings.showIntermediateResults ) showWithBdv( registeredImages, "registered" );

					Utils.log( "Creating projections..." );
					final ArrayList< ImagePlus > projections = createProjections( registeredImages );

					Utils.log( "Saving projections..." );
					saveImages( inputPath, projections );

					// Save full registered stack
					final RandomAccessibleInterval< T > transformedWithImagePlusDimensionOrder = Utils.copyAsArrayImg( Views.permute( registeredImages, 2, 3 ) );
					final ImagePlus transformedImagePlus = ImageJFunctions.wrap( transformedWithImagePlusDimensionOrder, "transformed" );
					final String outputPath = inputPath + "-registered.tif";
					Utils.log( "Saving registered image: " + outputPath );
					new FileSaver( transformedImagePlus ).saveAsTiff( outputPath );

					// Save svb non-registered projection
					RandomAccessibleInterval< T > shavenbaby = getShavenBabyImage( getImages( inputImagePlus ) );
					RandomAccessibleInterval shavenbabyMaximum = new Projection( shavenbaby, Z ).maximum();
					new FileSaver( ImageJFunctions.wrap( shavenbabyMaximum, "" ) ).saveAsTiff( inputPath + "-projection-ch1-raw.tif" );

					// Save ch2 non-registered projection
					RandomAccessibleInterval< T > ch2 = getChannel2Image( getImages( inputImagePlus ) );
					RandomAccessibleInterval ch2Maximum = new Projection( ch2, Z ).maximum();
					new FileSaver( ImageJFunctions.wrap( ch2Maximum, "" ) ).saveAsTiff( inputPath + "-projection-ch2-raw.tif" );

				}
			}
		}

		Utils.log( "Done!" );


	}

	public boolean acceptFile( String fileNameEndsWith, String file )
	{
		final String[] fileNameEndsWithList = fileNameEndsWith.split( "," );

		for ( String endsWith : fileNameEndsWithList )
		{
			if ( file.endsWith( endsWith.trim() ) )
			{
				return true;
			}
		}

		return false;
	}


	public void saveImages( String inputPath, ArrayList< ImagePlus > imps )
	{
		for ( ImagePlus imp : imps )
		{
			final String outputPath = inputPath + "-" + imp.getTitle() + ".tif";
			FileSaver fileSaver = new FileSaver( imp );
			fileSaver.saveAsTiff( outputPath );
		}
	}

	public ArrayList< ImagePlus > createProjections( RandomAccessibleInterval< T > images )
	{
		int Z = 2;

		ArrayList< ImagePlus > projections = new ArrayList<>(  );

		for ( int channelId = 0; channelId < images.dimension( 3 ); ++channelId )
		{
			RandomAccessibleInterval channel = Views.hyperSlice( images, 3, channelId );

			// top
			long rangeMin = (long) ( settings.finalProjectionMinDistanceToCenter / settings.outputResolution );
			long rangeMax = images.max( Z );
			Projection projection = new Projection( channel, Z, rangeMin, rangeMax );
			RandomAccessibleInterval maximum = projection.maximum();
			ImagePlus wrap = ImageJFunctions.wrap( maximum, "top-projection-ch" + ( channelId + 1 ) );
			projections.add( wrap );

			// bottom
			rangeMin = images.min( Z );
			rangeMax = - (long) ( settings.finalProjectionMinDistanceToCenter / settings.outputResolution );
			projection = new Projection( channel, Z, rangeMin, rangeMax );
			maximum = projection.maximum();
			wrap = ImageJFunctions.wrap( maximum, "bottom-projection-ch" + ( channelId + 1 ) );
			projections.add( wrap );

			// full
			rangeMin = images.min( Z );
			rangeMax = images.max( Z );
			projection = new Projection( channel, Z, rangeMin, rangeMax );
			maximum = projection.maximum();
			wrap = ImageJFunctions.wrap( maximum, "projection-ch" + ( channelId + 1 ) );
			projections.add( wrap );



		}

		return projections;
	}

	public void showWithBdv( RandomAccessibleInterval< T > transformed, String title )
	{
		Bdv bdv = BdvFunctions.show( transformed, title, BdvOptions.options().axisOrder( AxisOrder.XYZC ) );
		final ArrayList< RealPoint > points = new ArrayList<>();
		points.add( new RealPoint( new double[]{0,0,0} ));
		BdvFunctions.showPoints( points, "origin", BdvOptions.options().addTo( bdv ) );
	}

	public RandomAccessibleInterval< T > alignAndMaskImages(
			ImagePlus imagePlus,
			ShavenBabyRegistration registration )
	{
		final double[] inputCalibration = Utils.getCalibration( imagePlus );
		RandomAccessibleInterval< T > images = getImages( imagePlus );
		RandomAccessibleInterval< T > shavenbaby = getShavenBabyImage( images );
		RandomAccessibleInterval< T > channel2 = getChannel2Image( images );

		/**
		 * Compute registration
		 */
		Utils.log( "Computing registration...." );
		registration.run( shavenbaby, channel2, inputCalibration );

		/**
		 * Get transformation for demanded output resolution
		 */
		final AffineTransform3D registrationTransform = registration.getRegistrationTransform( inputCalibration, settings.outputResolution );
		if ( registrationTransform == null ) return null;

		/**
		 * Get transformation for demanded output resolution
		 */
		Utils.log( "Applying intensity correction to all channels...." );
		final RandomAccessibleInterval< T > correctedImages =
				createIntensityCorrectedImages(
						images,
						registration.getCorrectedCalibration()[ Z ],
						registration.getCoverslipPosition()  );

		Utils.log( "Creating registered and masked images (can take some time)..." );
		ArrayList< RandomAccessibleInterval< T > > registeredImages =
				Transforms.transformAllChannels(
						correctedImages,
						registrationTransform,
						settings.getOutputImageInterval()
				);

		// apply masking
		final RandomAccessibleInterval< BitType > alignedMaskAtOutputResolution = registration.createAlignedMask( settings.outputResolution, settings.getOutputImageInterval();
		registeredImages = Utils.maskAllChannels( registeredImages, alignedMaskAtOutputResolution );

		return Views.stack( registeredImages );
	}

	public RandomAccessibleInterval< T > createIntensityCorrectedImages( RandomAccessibleInterval< T > images,
																		 double axialCalibration,
																		 double coverslipPosition )
	{

		final RefractiveIndexMismatchCorrectionSettings correctionSettings = new RefractiveIndexMismatchCorrectionSettings();
		correctionSettings.pixelCalibrationMicrometer = axialCalibration;
		correctionSettings.coverslipPositionMicrometer = coverslipPosition;
		correctionSettings.intensityDecayLengthMicrometer = settings.refractiveIndexIntensityCorrectionDecayLength;

		return RefractiveIndexMismatchCorrections.createIntensityCorrectedImages( images, correctionSettings  );
	}


	public RandomAccessibleInterval< T > getImages( ImagePlus imagePlus )
	{
		RandomAccessibleInterval< T > images = ImageJFunctions.wrap( imagePlus );

		int numChannels = imagePlus.getNChannels();

		if ( numChannels == 1 )
		{
			Views.addDimension( images );
		}
		else
		{
			images = Views.permute( images, Utils.imagePlusChannelDimension, 3 );
		}

		return images;
	}

	private RandomAccessibleInterval< T > getShavenBabyImage( RandomAccessibleInterval< T > images )
	{
		RandomAccessibleInterval< T > rai = Views.hyperSlice( images, 3, settings.shavenbabyChannelIndexOneBased - 1 );

		return rai;
	}

	private RandomAccessibleInterval<T> getChannel2Image( RandomAccessibleInterval<T> images )
	{
		RandomAccessibleInterval< T > rai = Views.hyperSlice( images, 3, settings.amnioserosaChannelIndexOneBased - 1 );

		return rai;
	}

	public void setSettingsFromUI()
	{
		settings.showIntermediateResults = showIntermediateResults;
		settings.registrationResolution = registrationResolution;
		settings.closingRadius = 0;
		settings.outputResolution = outputResolution;
		//settings.refractiveIndexAxialCalibrationCorrectionFactor = refractiveIndexAxialCalibrationCorrectionFactor;
		settings.refractiveIndexIntensityCorrectionDecayLength = refractiveIndexIntensityCorrectionDecayLength;
		settings.thresholdModality = "";
		//settings.thresholdInUnitsOfBackgroundPeakHalfWidth = thresholdInUnitsOfBackgroundPeakHalfWidth;
		settings.rollAngleComputationMethod = rollAngleComputationMethod;
	}


}
