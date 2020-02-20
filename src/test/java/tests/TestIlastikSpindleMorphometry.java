package tests;

import de.embl.cba.morphometry.commands.SpindleMorphometryCommand;
import de.embl.cba.morphometry.spindle.SpindleMeasurements;
import de.embl.cba.morphometry.spindle.SpindleMorphometry;
import de.embl.cba.morphometry.spindle.SpindleMorphometrySettings;
import loci.common.DebugTools;
import net.imagej.ImageJ;
import net.imglib2.type.numeric.RealType;
import org.ilastik.ilastik4ij.ui.IlastikOptions;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;

public class TestIlastikSpindleMorphometry
{
	private static boolean showOutput = false;
	private static boolean showIntermediateResults = false;

	//	@Test
	public < R extends RealType< R > > void run()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final IlastikOptions ilastikOptions = new IlastikOptions();
		ilastikOptions.setExecutableFile(new File("/Applications/ilastik-1.3.3-OSX.app/Contents/MacOS/ilastik"));

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
		command.opService = ij.op();
		command.datasetService = ij.dataset();
		command.statusService = ij.status();
		command.logService = ij.log();
		command.optionsService = ij.options();

		command.inputImageFile = new File(
				TestIlastikSpindleMorphometry.class.getResource(
						"../test-data/spindle/SpindleWidthSmall.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.voxelSpacingDuringAnalysis = 0.25;
		command.showIntermediateResults = false;
		command.ilastikOptions = ilastikOptions;
		command.classifier = SpindleMorphometry.ILASTIK;
		command.classifierFile = new File("/Users/tischer/Documents/tobias-kletter/2020-02-ilastik-test/20191206_DNA_Segmentation_2Ch.ilp");
		command.saveResults = true;
		command.settings.showOutputImage = showOutput;
		command.outputDirectory = new File("/Users/tischer/" +
				"Documents/fiji-plugin-morphometry/src" +
				"/test/resources/test-data/spindle/output" );
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidth = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		final Double dnaLateralExtend = ( Double ) measurements.get( 0 ).get(
				SpindleMeasurements.getDnaLateralExtendKey() );

		//assertEquals( 9.5, dnaLateralExtend, 1.0 );
		//assertEquals( 9.0, spindleWidth, 1.0 );
	}

	public static void main( String[] args )
	{
		showOutput = false;
		showIntermediateResults = false;

		final ImageJ imageJ = new ImageJ();
		imageJ.ui().showUI();

		new TestIlastikSpindleMorphometry().run();
	}

}
