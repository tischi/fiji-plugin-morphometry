package tests;

import de.embl.cba.morphometry.commands.SpindleMorphometryCommand;
import de.embl.cba.morphometry.spindle.SpindleMeasurements;
import loci.common.DebugTools;
import net.imagej.ImageJ;
import net.imglib2.type.numeric.RealType;
import org.junit.Test;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;

public class TestSpindleMorphometry
{

	public static boolean showImageJUI = false;

	@Test
	public < R extends RealType< R > > void testSmallSpindle( )
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		if ( showImageJUI )
			ij.ui().showUI();

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/SpindleWidthSmall.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.voxelSpacingDuringAnalysis = 0.25;
		command.showIntermediateResults = false;
		command.saveResults = true;
		command.settings.showOutputImage = true;
		command.outputDirectory = new File("/Users/tischer/" +
				"Documents/fiji-plugin-morphometry/src" +
				"/test/resources/test-data/spindle/output" );
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidth = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		final Double dnaLateralExtend = ( Double ) measurements.get( 0 ).get(
				getDnaLateralExtendKey() );

		assertEquals( 9.5, dnaLateralExtend, 1.0 );
		assertEquals( 9.0, spindleWidth, 1.0 );
	}

	@Test
	public < R extends RealType< R > > void testLargeSpindle()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/SpindleWidthLarge.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.showIntermediateResults = false;
		command.saveResults = false;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidth = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		final Double dnaLateralExtend = ( Double ) measurements.get( 0 ).get(
				getDnaLateralExtendKey() );

		assertEquals( 12.0, dnaLateralExtend, 1.0 );
		assertEquals( 11.0, spindleWidth, 1.0 );
	}


	@Test
	public < R extends RealType< R > > void testSpindleWithBrightOtherDNA()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/BrightOtherDNA.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 2;
		command.dnaChannelIndexOneBased = 1;
		command.showIntermediateResults = false;
		command.saveResults = false;
		command.outputDirectory = new File("/Users/tischer/" +
				"Documents/fiji-plugin-morphometry/src" +
				"/test/resources/test-data/spindle/output" );
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidth = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		final Double dnaLateralExtend = ( Double ) measurements.get( 0 ).get(
				getDnaLateralExtendKey() );

		assertEquals( 15.5, dnaLateralExtend, 1.0 );
		assertEquals( 13.0, spindleWidth, 1.0 );
	}

	private String getDnaLateralExtendKey()
	{
		return SpindleMeasurements.DNA_LATERAL_EXTEND
				+ SpindleMeasurements.SEP + SpindleMeasurements.LENGTH_UNIT;
	}

	@Test
	public < R extends RealType< R > > void testDimDNA()
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();
		//ij.ui().showUI();

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
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


	public static void main( String[] args )
	{
//		new TestSpindleMorphometry().testDimDNA();
//		new TestSpindleMorphometry().testSpindleWithBrightOtherDNA();
		showImageJUI = true; new TestSpindleMorphometry().testSmallSpindle();
//		new TestSpindleMorphometry().testLargeSpindle();
	}

}
