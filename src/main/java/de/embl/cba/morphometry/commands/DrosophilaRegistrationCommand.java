package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Projection;
import de.embl.cba.morphometry.drosophila.registration.DrosophilaRegistrationSettings;
import de.embl.cba.morphometry.drosophila.registration.DrosophilaSingleChannelRegistration;
import de.embl.cba.morphometry.refractiveindexmismatch.RefractiveIndexMismatchCorrectionSettings;
import de.embl.cba.morphometry.refractiveindexmismatch.RefractiveIndexMismatchCorrections;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.transforms.utils.Transforms;
import ij.ImagePlus;
import ij.io.FileSaver;
import net.imagej.DatasetService;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
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


@Plugin(type = Command.class, menuPath = "Plugins>Registration>EMBL>Drosophila Registration" )
public class DrosophilaRegistrationCommand < T extends RealType< T > & NativeType< T > > implements Command
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

//	@Parameter( choices = { DrosophilaRegistrationSettings.MOMENTS })
	public String longAxisAngleAlignmentMethod = DrosophilaRegistrationSettings.MOMENTS;

//	@Parameter( choices = { DrosophilaRegistrationSettings.INTENSITY })
	public String longAxisFlippingAlignmentMethod = DrosophilaRegistrationSettings.INTENSITY;

//	@Parameter( choices = { DrosophilaRegistrationSettings.INTENSITY })
	public String rollAngleAlignmentMethod = DrosophilaRegistrationSettings.INTENSITY;

	@Parameter
	public int alignmentChannelIndexOneBased = settings.alignmentChannelIndexOneBased;

//	@Parameter
//	public int secondaryChannelIndexOneBased = settings.secondaryChannelIndexOneBased;

	@Parameter
	public double registrationResolution = settings.registrationResolution;

	@Parameter
	public double outputResolution = settings.outputResolution;

	@Parameter
	public double refractiveIndexIntensityCorrectionDecayLength = settings.refractiveIndexIntensityCorrectionDecayLength;


	public void run()
	{
		setSettingsFromUI();

		final DrosophilaSingleChannelRegistration registration =
				new DrosophilaSingleChannelRegistration( settings, opService );

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
					Logger.log( " " );
					Logger.log( "Reading: " + inputPath + "..." );
					final ImagePlus inputImagePlus = openWithBioFormats( inputPath );

					if ( inputImagePlus == null )
					{
						logService.error( "Error opening labelMaskFile: " + inputPath );
						continue;
					}

					/**
					 * Register
					 */

					RandomAccessibleInterval< T > registeredImages =
							createAlignedImages( inputImagePlus, registration );

					if ( registeredImages == null )
					{
						Logger.log( "ERROR: Could not find central embryo" );
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
//					RandomAccessibleInterval< T > channel1Image = getChannelImage( getChannelImages( inputImagePlus ) );
//					RandomAccessibleInterval shavenbabyMaximum = new Projection( channel1Image, Z ).maximum();
//					new FileSaver( ImageJFunctions.wrap( shavenbabyMaximum, "" ) ).saveAsTiff( outputFilePathStump + "-projection-ch1-raw.tif" );
//
//					// Save ch2 non-registered projection
//					RandomAccessibleInterval< T > channel2Image = getChannel2Image( getChannelImages( inputImagePlus ) );
//					RandomAccessibleInterval ch2Maximum = new Projection( channel2Image, Z ).maximum();
//					new FileSaver( ImageJFunctions.wrap( ch2Maximum, "" ) ).saveAsTiff( outputFilePathStump + "-projection-ch2-raw.tif" );

				}
			}
		}

		Logger.log( "Done!" );


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
		Logger.log( "Saving registered image: " + outputPath );
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
			RandomAccessibleInterval channel = ImageIO.getChannelImage( images, channelId );

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

	public RandomAccessibleInterval< T > createAlignedImages(
			ImagePlus imagePlus,
			DrosophilaSingleChannelRegistration registration )
	{
		final double[] inputCalibration = Utils.getCalibration( imagePlus );
		RandomAccessibleInterval< T > images = ImageIO.getChannelImages( imagePlus );
		RandomAccessibleInterval< T > image = ImageIO.getChannelImage( images, alignmentChannelIndexOneBased - 1  );

		/**
		 * Compute registration
		 */
		Logger.log( "Computing registration...." );
		registration.run( image, inputCalibration );

		/**
		 * Apply intensity correction
		 */
		Logger.log( "Applying intensity correction to all channels...." );
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
		Logger.log( "Creating registered and masked images (can take some time)..." );
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
				= registration.getAlignedMask( settings.outputResolution, settings.getOutputImageInterval() );
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


	public void setSettingsFromUI()
	{
		settings.showIntermediateResults = showIntermediateResults;
		settings.registrationResolution = registrationResolution;
		settings.outputResolution = outputResolution;
		settings.refractiveIndexIntensityCorrectionDecayLength = refractiveIndexIntensityCorrectionDecayLength;
		settings.thresholdModality = "";
		settings.rollAngleComputationMethod = rollAngleAlignmentMethod;
		settings.alignmentChannelIndexOneBased = alignmentChannelIndexOneBased;
	}


}
