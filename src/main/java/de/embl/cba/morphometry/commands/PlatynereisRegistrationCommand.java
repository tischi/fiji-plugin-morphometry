package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.registration.platynereis.PlatynereisRegistration;
import de.embl.cba.morphometry.registration.platynereis.PlatynereisRegistrationSettings;
import de.embl.cba.transforms.utils.Transforms;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.plugin.FolderOpener;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import java.io.File;
import java.util.ArrayList;


@Plugin(type = Command.class, menuPath = "Plugins>Registration>EMBL>Platynereis Registration" )
public class PlatynereisRegistrationCommand< R extends RealType< R > & NativeType< R > > implements Command
{
	public static final String SAVE_AS_BDV = "Save as bdv";

	@Parameter
	public OpService opService;

	PlatynereisRegistrationSettings settings = new PlatynereisRegistrationSettings();

	@Parameter( style = "directory" )
	public File inputDirectory;

	// TODO: change to path
	@Parameter( style = "directory" )
	public File outputDirectory;

	@Parameter( label = "Only open files ending with" )
	public String fileNameEndsWith = ".tiff";

	@Parameter( label = "Show intermediate results (for debugging)")
	public boolean showIntermediateResults = settings.showIntermediateResults;

	@Parameter( label = "Input image isotropic pixel size [micrometer]" )
	public double inputResolutionMicrometer;

	@Parameter( label = "Invert image (pixels of interest are dark)" )
	public boolean invertImage = settings.invertImage;

	@Parameter( label = "Registration resolution [micrometer]")
	public double registrationResolution = settings.registrationResolution;

	@Parameter( label = "Output resolution [micrometer]" )
	public double outputResolution = settings.outputResolution;

	@Parameter( label = "Output modality", choices = { SAVE_AS_BDV } )
	public String outputModality = SAVE_AS_BDV;


	public void run()
	{
		setSettings();

		// Open image
		//
		Logger.log( "Opening image: " + inputDirectory.getAbsolutePath() );
		final ImagePlus imagePlus = FolderOpener.open(
				inputDirectory.getAbsolutePath(),
				" file=(.*" + fileNameEndsWith + ")" );
		final RandomAccessibleInterval< R > channelImages = Utils.getChannelImages( imagePlus );
		RandomAccessibleInterval< R > image = Utils.getChannelImage( channelImages, 0 );

		// Find registration
		//
		Logger.log( "Registering..." );
		final PlatynereisRegistration< R > registration = new PlatynereisRegistration<>( settings, opService );
		registration.run( image );

		// Apply registration
		//
		Logger.log( "Preparing registered images..." );
		ArrayList< RandomAccessibleInterval< R > > registeredImages =
				Transforms.transformAllChannels(
						channelImages,
						registration.getRegistrationTransform( settings.inputCalibration, settings.outputResolution ),
						settings.getOutputImageInterval()
				);


		// Save registered image
		//
		final String outputFilePathStump = outputDirectory + File.separator + inputDirectory.getName();

		saveRegisteredImage( outputFilePathStump, registeredImages.get( 0 ) );

		Logger.log( "Done!" );
	}

	public void saveRegisteredImage( String outputFilePathStump, RandomAccessibleInterval< R > registeredImage )
	{
		registeredImage = Utils.get3DRaiAs5DRaiWithImagePlusDimensionOrder( registeredImage );

		final ImagePlus registered = ImageJFunctions.wrap( registeredImage, "transformed" );
		registered.getCalibration().setUnit( "micrometer" );
		registered.getCalibration().pixelWidth = settings.outputResolution;
		registered.getCalibration().pixelHeight = settings.outputResolution;
		registered.getCalibration().pixelDepth = settings.outputResolution;

		final String outputPath = outputFilePathStump + "-aligned.tif";
		Logger.log( "Applying registration and saving registered image: " + outputPath );
		new FileSaver( registered ).saveAsTiff( outputPath );
	}


	public void setSettings()
	{
		settings.showIntermediateResults = showIntermediateResults;
		settings.registrationResolution = registrationResolution;
		settings.outputResolution = outputResolution;
		settings.invertImage = invertImage;
		settings.inputCalibration = new double[]{
				inputResolutionMicrometer,
				inputResolutionMicrometer,
				inputResolutionMicrometer };

	}


}
