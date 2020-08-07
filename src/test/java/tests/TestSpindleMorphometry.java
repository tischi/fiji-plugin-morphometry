package tests;

import de.embl.cba.morphometry.commands.Spindle3DCommand;
import de.embl.cba.morphometry.spindle.SpindleMeasurements;
import de.embl.cba.morphometry.spindle.SpindleMorphometrySettings;
import loci.common.DebugTools;
import net.imagej.ImageJ;
import net.imglib2.type.numeric.RealType;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;

public class TestSpindleMorphometry
{
	private static boolean showOutput = false;
	private static boolean showIntermediateResults = false;


	//	@Test
	public < R extends RealType< R > > void testSmallSpindle( )
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final Spindle3DCommand< R > command = new Spindle3DCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/SpindleWidthSmall.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.voxelSpacingDuringAnalysis = 0.25;
		command.showIntermediateResults = false;
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

		assertEquals( 9.5, dnaLateralExtend, 1.0 );
		assertEquals( 9.0, spindleWidth, 1.0 );
	}

//	@Test
	public < R extends RealType< R > > void testExceptionCausingImage( )
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final Spindle3DCommand< R > command = new Spindle3DCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/ExceptionTest00.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.voxelSpacingDuringAnalysis = 0.25;
		command.showIntermediateResults = false;
		command.saveResults = false;
		command.settings.showOutputImage = showOutput;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidth = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		final Double dnaLateralExtend = ( Double ) measurements.get( 0 ).get(
				SpindleMeasurements.getDnaLateralExtendKey() );

//		assertEquals( 9.5, dnaLateralExtend, 1.0 );
//		assertEquals( 9.0, spindleWidth, 1.0 );
	}

//	@Test
	public < R extends RealType< R > > void testLargeSpindle()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final Spindle3DCommand< R > command = new Spindle3DCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/SpindleWidthLarge.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.showIntermediateResults = false;
		command.saveResults = false;
		command.cellCenterDetectionMethodChoice = SpindleMorphometrySettings.CellCenterDetectionMethod.None.toString();
		command.settings.showOutputImage = showOutput;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidthMax = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		final Double dnaLateralExtend = ( Double ) measurements.get( 0 ).get(
				SpindleMeasurements.getDnaLateralExtendKey() );

		final Double spindleVolume = ( Double ) measurements.get( 0 ).get( SpindleMeasurements.getSpindleVolumeKey() );

		assertEquals( 12.0, dnaLateralExtend, 1.0 );
		assertEquals( 11.0, spindleWidthMax, 1.0 );
		assertEquals( 460, spindleVolume, 50 );
	}

	//@Test
	public < R extends RealType< R > > void testSpindleVolume00()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final Spindle3DCommand< R > command = new Spindle3DCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/SpindleVolumeTest00.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.showIntermediateResults = false;
		command.saveResults = false;
		command.settings.showOutputImage = showOutput;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidth = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		final Double dnaLateralExtend = ( Double ) measurements.get( 0 ).get(
				SpindleMeasurements.getDnaLateralExtendKey() );

		final Double spindleVolume = ( Double ) measurements.get( 0 ).get( SpindleMeasurements.getSpindleVolumeKey() );

//		assertEquals( 12.0, dnaLateralExtend, 1.0 );
//		assertEquals( 11.0, spindleWidth, 1.0 );
		assertEquals( 330, spindleVolume, 50 );
	}

	//@Test
	public < R extends RealType< R > > void testSpindleVolume01()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final Spindle3DCommand< R > command = new Spindle3DCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/SpindleVolumeTest01.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.showIntermediateResults = false;
		command.saveResults = false;
		command.settings.showOutputImage = showOutput;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidth = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		final Double dnaLateralExtend = ( Double ) measurements.get( 0 ).get(
				SpindleMeasurements.getDnaLateralExtendKey() );

		final Double spindleVolume = ( Double ) measurements.get( 0 ).get( SpindleMeasurements.getSpindleVolumeKey() );

//		assertEquals( 12.0, dnaLateralExtend, 1.0 );
//		assertEquals( 11.0, spindleWidth, 1.0 );
		assertEquals( 400, spindleVolume, 50 );
	}

	//@Test
	public < R extends RealType< R > > void testSpindleVolume02()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final Spindle3DCommand< R > command = new Spindle3DCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/SpindleVolumeTest02.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.showIntermediateResults = false;
		command.saveResults = false;
		command.settings.showOutputImage = showOutput;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidth = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		final Double dnaLateralExtend = ( Double ) measurements.get( 0 ).get(
				SpindleMeasurements.getDnaLateralExtendKey() );

		final Double spindleVolume = ( Double ) measurements.get( 0 ).get( SpindleMeasurements.getSpindleVolumeKey() );

//		assertEquals( 12.0, dnaLateralExtend, 1.0 );
//		assertEquals( 11.0, spindleWidth, 1.0 );
		assertEquals( 400, spindleVolume, 50 );
	}

//	@Test
	public < R extends RealType< R > > void weirdDNA3Channels()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final Spindle3DCommand< R > command = new Spindle3DCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/WeirdDNA3Channels.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 2;
		command.dnaChannelIndexOneBased = 1;
		command.showIntermediateResults = false;
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

//		assertEquals( 12.0, dnaLateralExtend, 1.0 );
//		assertEquals( 11.0, spindleWidth, 1.0 );
	}

//	@Test
	public < R extends RealType< R > > void testSpindleWithBrightOtherDNA()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final Spindle3DCommand< R > command = new Spindle3DCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/BrightOtherDNA.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 2;
		command.dnaChannelIndexOneBased = 1;
		command.cellCenterDetectionMethodChoice = SpindleMorphometrySettings.CCDM_NONE;
		command.showIntermediateResults = showIntermediateResults;
		command.saveResults = false;
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

		assertEquals( 15.5, dnaLateralExtend, 1.0 );
		assertEquals( 13.0, spindleWidth, 1.0 );
	}

	//@Test
	public < R extends RealType< R > > void testDimDNA()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();
		//ij.ui().showUI();

		final Spindle3DCommand< R > command = new Spindle3DCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/DimDNA.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.showIntermediateResults = false;
		command.saveResults = false;
		command.outputDirectory = new File("/Users/tischer/" +
				"Documents/fiji-plugin-morphometry/src" +
				"/test/resources/test-data/spindle/output" );
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final String comment = ( String) measurements.get( 0 ).get(
				SpindleMeasurements.COMMENT );

		assertEquals( comment, SpindleMeasurements.ANALYSIS_INTERRUPTED_LOW_DYNAMIC_DNA );
	}

	//@Test
	public < R extends RealType< R > > void testLargeCoV()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();
		//ij.ui().showUI();

		final Spindle3DCommand< R > command = new Spindle3DCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/HeterogeneousSpindle.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 3;
		command.dnaChannelIndexOneBased = 1;
		command.showIntermediateResults = false;
		command.settings.showOutputImage = showOutput;
		command.saveResults = false;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final double spindleCoV = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.SPINDLE_TUBULIN_COV );

		assertEquals( 1.22, spindleCoV, 0.2 );
	}

	//@Test
	public < R extends RealType< R > > void testSmallCoV()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();
		//ij.ui().showUI();

		final Spindle3DCommand< R > command = new Spindle3DCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/HomogenousSpindle.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.showIntermediateResults = false;
		command.settings.showOutputImage = showOutput;
		command.saveResults = false;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final double spindleCoV = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.SPINDLE_TUBULIN_COV );

		final double spindleLength = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleLengthKey() );

		assertEquals( 0.73, spindleCoV, 0.2 );
		assertEquals( 10.0, spindleLength, 1.0 );
	}

//	@Test
	public < R extends RealType< R > > void testSpindleLength()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();
		//ij.ui().showUI();

		final Spindle3DCommand< R > command = new Spindle3DCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/SpindleLength9.5.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 3;
		command.dnaChannelIndexOneBased = 1;
		command.showIntermediateResults = false;
		command.settings.showOutputImage = showOutput;
		command.saveResults = false;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final double spindleLength = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleLengthKey() );

		// manual GT is 10.8
		assertEquals( 10.8, spindleLength, 1.0 );
	}

	public static void main( String[] args )
	{
		showOutput = false;
		showIntermediateResults = false;

		final ImageJ imageJ = new ImageJ();
		imageJ.ui().showUI();

		new TestSpindleMorphometry().testLargeSpindle();
//		new TestSpindleMorphometry().testDimDNA();
//		new TestSpindleMorphometry().testSpindleWithBrightOtherDNA();
//		showImageJUI = true; new TestSpindleMorphometry().testSmallSpindle();
//		showImageJUI = true; new TestSpindleMorphometry().weirdDNA3Channels();
//		new TestSpindleMorphometry().testSmallCoV();
//		new TestSpindleMorphometry().testSpindleLength();
	}

}
