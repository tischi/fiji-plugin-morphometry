package tests;

import de.embl.cba.morphometry.commands.SpindleMorphometryCommand;
import de.embl.cba.morphometry.spindle.SpindleMeasurements;
import net.imagej.ImageJ;
import net.imglib2.type.numeric.RealType;
import org.junit.Test;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;

public class TestSpindleMorphometry
{
	@Test
	public < R extends RealType< R > > void testSmallSpindle( )
	{
		final ImageJ ij = new ImageJ();

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/SpindleWidthSmall.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.voxelSpacingDuringAnalysis = 0.25;
		command.showIntermediateResults = false;
		command.saveResults = false;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidth = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		assertEquals( spindleWidth, 7.5, 1.0 );
	}

	@Test
	public < R extends RealType< R > > void testLargeSpindle()
	{
		final ImageJ ij = new ImageJ();

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/SpindleWidthLarge.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.voxelSpacingDuringAnalysis = 0.25;
		command.showIntermediateResults = false;
		command.saveResults = false;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidth = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		assertEquals( spindleWidth, 8.79, 1.0 );
	}


	@Test
	public < R extends RealType< R > > void testSpindleWithBrightOtherDNA()
	{
		final ImageJ ij = new ImageJ();

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/BrightOtherDNA.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 2; // This is the other way around than in the others!
		command.dnaChannelIndexOneBased = 1;
		command.voxelSpacingDuringAnalysis = 0.25;
		command.showIntermediateResults = false;
		command.saveResults = false;
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final Double spindleWidth = ( Double) measurements.get( 0 ).get(
				SpindleMeasurements.getSpindleWidthMaxKey() );

		assertEquals( spindleWidth, 9.25, 1.0 );
	}


	@Test
	public < R extends RealType< R > > void testDimDNA()
	{
		final ImageJ ij = new ImageJ();
		//ij.ui().showUI();

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File(
				TestSpindleMorphometry.class.getResource(
						"../test-data/spindle/DimDNA.zip" ).getFile() );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;
		command.voxelSpacingDuringAnalysis = 0.25;
		command.minimalDynamicRange = 20;
		command.showIntermediateResults = false;
		command.saveResults = true;
		command.outputDirectory = new File("/Users/tischer/" +
				"Documents/fiji-plugin-morphometry/src" +
				"/test/resources/test-data/spindle/output" );
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();

		final String comment = ( String) measurements.get( 0 ).get(
				SpindleMeasurements.COMMENT );

		assertEquals( comment, SpindleMeasurements.TOO_LOW_DYNAMIC_RANGE_IN_DNA_IMAGE  );
	}


	public static void main( String[] args )
	{
		new TestSpindleMorphometry().testDimDNA();
//		new TestSpindleMorphometry().testSpindleWithBrightOtherDNA();
//		new TestSpindleMorphometry().testSmallSpindle();
//		new TestSpindleMorphometry().testLargeSpindle();
	}

}
