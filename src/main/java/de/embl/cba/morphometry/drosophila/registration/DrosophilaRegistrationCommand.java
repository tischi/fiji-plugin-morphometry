package de.embl.cba.morphometry.drosophila.registration;

import bdv.util.*;
import de.embl.cba.morphometry.Projection;
import de.embl.cba.morphometry.refractiveindexmismatch.RefractiveIndexMismatchCorrectionSettings;
import de.embl.cba.morphometry.refractiveindexmismatch.RefractiveIndexMismatchCorrections;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.transforms.utils.Transforms;
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
public class DrosophilaRegistrationCommand<T extends RealType<T> & NativeType< T > > implements Command
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

	DrosophilaRegistrationSettings settings = new DrosophilaRegistrationSettings();

	public static final String FROM_DIRECTORY = "From directory";
	public static final String CURRENT_IMAGE = "Current image";

	@Parameter( choices = { FROM_DIRECTORY, CURRENT_IMAGE })
	public String inputModality = FROM_DIRECTORY;

	@Parameter( style = "directory" )
	public File inputDirectory;

	@Parameter( style = "directory" )
	public File outputDirectory;

	@Parameter
	public String fileNameEndsWith = ".czi,.lsm";

	@Parameter
	public boolean showIntermediateResults = settings.showIntermediateResults;

	@Parameter( choices = { DrosophilaRegistrationSettings.PRIMARY_CHANNEL_MOMENTS })
	public String longAxisAngleAlignmentMethod = DrosophilaRegistrationSettings.PRIMARY_CHANNEL_MOMENTS;

	@Parameter( choices = { DrosophilaRegistrationSettings.PRIMARY_CHANNEL_INTENSITY })
	public String longAxisFlippingAlignmentMethod = DrosophilaRegistrationSettings.PRIMARY_CHANNEL_INTENSITY;

	@Parameter( choices = {
			DrosophilaRegistrationSettings.CENTROID_SHAPE_BASED_ROLL_TRANSFORM,
			DrosophilaRegistrationSettings.PROJECTION_SHAPE_BASED_ROLL_TRANSFORM,
			DrosophilaRegistrationSettings.SECONDARY_CHANNEL_INTENSITY })
	public String rollAngleAlignmentMethod = DrosophilaRegistrationSettings.SECONDARY_CHANNEL_INTENSITY;


	@Parameter
	public int primaryChannelIndexOneBased = settings.primaryChannelIndexOneBased;

	@Parameter
	public int secondaryChannelIndexOneBased = settings.secondaryChannelIndexOneBased;

	@Parameter
	public double registrationResolution = settings.registrationResolution;

	@Parameter
	public double outputResolution = settings.outputResolution;

	@Parameter
	public double refractiveIndexIntensityCorrectionDecayLength = settings.refractiveIndexIntensityCorrectionDecayLength;


	public void run()
	{
		setSettingsFromUI();

		final DrosphilaRegistration registration = new DrosphilaRegistration( settings, opService );

		if ( inputModality.equals( CURRENT_IMAGE ) && imagePlus != null )
		{
//			RandomAccessibleInterval< T > transformed = createAlignedImages( imagePlus, registration );
//			showWithBdv( transformed, "registered" );
//			ImageJFunctions.show( Views.permute( transformed, 2, 3 ) );
		}

		// TODO: move the batching out of this and use generic batching?
		if ( inputModality.equals( FROM_DIRECTORY ) )
		{
			String[] files = inputDirectory.list();

			for( String file : files )
			{
				if ( acceptFile( fileNameEndsWith, file ) )
				{
					final String outputFilePathStump = outputDirectory + File.separator + file;

					Utils.setNewLogFilePath( outputFilePathStump + ".log.txt" );

					/**
					 * Open the images
					 */

					final String inputPath = inputDirectory + File.separator + file;
					Utils.log( " " );
					Utils.log( "Reading: " + inputPath + "..." );
					final ImagePlus inputImagePlus = openWithBioFormats( inputPath );

					if ( inputImagePlus == null )
					{
						logService.error( "Error opening inputLabelMaskFile: " + inputPath );
						continue;
					}

					/**
					 * Register
					 */

					RandomAccessibleInterval< T > registeredImages = createAlignedImages( inputImagePlus, registration );

					if ( registeredImages == null )
					{
						Utils.log( "ERROR: Could not find central embryo" );
						continue;
					}

					/**
					 * Save registered images
					 */

					saveResults( outputFilePathStump, registeredImages );


					// other stuff
					//
//					RandomAccessibleInterval< T > watershed = (RandomAccessibleInterval) registration.getWatershedLabelImg();
//					new FileSaver( ImageJFunctions.wrap( watershed, "" ) ).saveAsTiff( outputFilePathStump + "-watershed.tif" );
//
//					Utils.log( "Creating projections..." );
//					final ArrayList< ImagePlus > projections = createProjections( registeredImages );
//
//					Utils.log( "Saving projections..." );
//					saveImages( outputFilePathStump, projections );
//

//					// Save ch1 non-registered projection
//					RandomAccessibleInterval< T > channel1Image = getChannel1Image( getImages( inputImagePlus ) );
//					RandomAccessibleInterval shavenbabyMaximum = new Projection( channel1Image, Z ).maximum();
//					new FileSaver( ImageJFunctions.wrap( shavenbabyMaximum, "" ) ).saveAsTiff( outputFilePathStump + "-projection-ch1-raw.tif" );
//
//					// Save ch2 non-registered projection
//					RandomAccessibleInterval< T > channel2Image = getChannel2Image( getImages( inputImagePlus ) );
//					RandomAccessibleInterval ch2Maximum = new Projection( channel2Image, Z ).maximum();
//					new FileSaver( ImageJFunctions.wrap( ch2Maximum, "" ) ).saveAsTiff( outputFilePathStump + "-projection-ch2-raw.tif" );

				}
			}
		}

		Utils.log( "Done!" );


	}

	public void saveResults( String outputFilePathStump, RandomAccessibleInterval< T > registeredImages )
	{
		// 3D image stack
		//
		final RandomAccessibleInterval< T > registeredWithImagePlusDimensionOrder = Utils.copyAsArrayImg( Views.permute( registeredImages, 2, 3 ) );
		final ImagePlus registered = ImageJFunctions.wrap( registeredWithImagePlusDimensionOrder, "transformed" );
		registered.getCalibration().setUnit( "micrometer" );
		registered.getCalibration().pixelWidth = settings.outputResolution;
		registered.getCalibration().pixelHeight = settings.outputResolution;
		registered.getCalibration().pixelDepth = settings.outputResolution;

		final String outputPath = outputFilePathStump + "-registered.tif";
		Utils.log( "Saving registered image: " + outputPath );
		new FileSaver( registered ).saveAsTiff( outputPath );
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


	public void saveImages( String outputPath, ArrayList< ImagePlus > imps )
	{
		for ( ImagePlus imp : imps )
		{
			final String outputPath2 = outputPath + "-" + imp.getTitle() + ".tif";
			FileSaver fileSaver = new FileSaver( imp );
			fileSaver.saveAsTiff( outputPath2 );
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

	public RandomAccessibleInterval< T > createAlignedImages(
			ImagePlus imagePlus,
			DrosphilaRegistration registration )
	{
		final double[] inputCalibration = Utils.getCalibration( imagePlus );
		RandomAccessibleInterval< T > images = getImages( imagePlus );
		RandomAccessibleInterval< T > channel1 = getChannel1Image( images );
		RandomAccessibleInterval< T > channel2 = getChannel2Image( images );

		/**
		 * Compute registration
		 */
		Utils.log( "Computing registration...." );
		registration.run( channel1, channel2, inputCalibration );

		/**
		 * Apply intensity correction
		 */
		Utils.log( "Applying intensity correction to all channels...." );
		final RandomAccessibleInterval< T > intensityCorrectedImages =
				createIntensityCorrectedImages(
						images,
						registration.getCorrectedCalibration()[ Z ],
						registration.getCoverslipPosition()  );


		/**
		 * Create transformation for desired output resolution
		 */
		final AffineTransform3D registrationTransform =
				registration.getRegistrationTransform(
						registration.getCorrectedCalibration(), settings.outputResolution );
		if ( registrationTransform == null ) return null;

		/**
		 * Apply transformation
		 */
		Utils.log( "Creating registered and masked images (can take some time)..." );
		ArrayList< RandomAccessibleInterval< T > > registeredImages =
				Transforms.transformAllChannels(
						intensityCorrectedImages,
						registrationTransform,
						settings.getOutputImageInterval()
				);

		/**
		 * Apply masking ( in order to remove other, potentially touching, embryos )
		 */
		final RandomAccessibleInterval< BitType > alignedMaskAtOutputResolution
				= registration.createAlignedMask( settings.outputResolution, settings.getOutputImageInterval() );
		registeredImages = Utils.maskAllChannels( registeredImages, alignedMaskAtOutputResolution, settings.showIntermediateResults );

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

	private RandomAccessibleInterval< T > getChannel1Image( RandomAccessibleInterval< T > images )
	{
		RandomAccessibleInterval< T > rai = Views.hyperSlice( images, 3, settings.primaryChannelIndexOneBased - 1 );

		return rai;
	}

	private RandomAccessibleInterval<T> getChannel2Image( RandomAccessibleInterval<T> images )
	{
		RandomAccessibleInterval< T > rai = Views.hyperSlice( images, 3, settings.secondaryChannelIndexOneBased - 1 );

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
		settings.rollAngleComputationMethod = rollAngleAlignmentMethod;
		settings.secondaryChannelIndexOneBased = secondaryChannelIndexOneBased;
		settings.primaryChannelIndexOneBased = primaryChannelIndexOneBased;
	}


}
