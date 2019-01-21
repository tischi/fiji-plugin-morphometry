package de.embl.cba.morphometry.microglia;

import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.skeleton.SkeletonCreator;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class MicrogliaMorphometry < T extends RealType< T > & NativeType< T > >
{

	private final ArrayList< RandomAccessibleInterval< T > > labelMaps;
	private final OpService opService;
	private ArrayList< HashMap< Integer, Map< String, Object > > > measurementsTimepointList;
	private ArrayList< RandomAccessibleInterval< BitType > > skeletons;

	public MicrogliaMorphometry( ArrayList< RandomAccessibleInterval< T > > labelMaps,
								 OpService opService )
	{
		this.labelMaps = labelMaps;
		this.opService = opService;
	}

	public void run()
	{
		createSkeletons( );

		measurementsTimepointList = Measurements.initMeasurements( labelMaps.size() );

		performMeasurements( );
	}

	private void performMeasurements( )
	{
		for ( int t = 0; t < labelMaps.size(); ++t )
		{
			final HashMap< Integer, Map< String, Object > > measurements = measurementsTimepointList.get( t );

			final ImgLabeling< Integer, IntType > imgLabeling = Utils.labelMapAsImgLabelingRobert( labelMaps.get( t ) );

			Measurements.measurePositions(
					measurements,
					imgLabeling,
					null);

			// Volumes ( = areas )
			Measurements.measureVolumes(
					measurements,
					imgLabeling);

			// Surfaces ( = perimeters )
			Measurements.measureSurface(
					measurements,
					imgLabeling,
					opService );

			// TODO: move to skeletonAnalyzer?
			Measurements.measureSkeletons(
					measurements,
					imgLabeling,
					skeletons.get( t ),
					opService );

			// Form factor could be calculated later, e.g. in R

			// Analyze Skeletons: length, branch-points, branches
			// avg branch-length = length / branches

			// Measure: distance travelled

			// Also,we are presently using MtrackJ to calculate velocity, distance travelled and displacement.
			// => I would recommend you do this in Excel as this is downstream analysis.
			// => What I can work on is a tool to upload your extended table again and view it on top of the objects

			// With the segmented microglia movie generated with your plugin, can we do automatic tracking?
			// The cells are already tracked.

			//	2-The next  challenge would be to measure phagocytosis of green particles and quantify "black holes" as we discussed last summer.


		}
	}

	private void createSkeletons( )
	{
		final SkeletonCreator skeletonCreator = new SkeletonCreator(
				Utils.labelMapsAsMasks( labelMaps ),
				opService );

		skeletonCreator.run();

		skeletons = skeletonCreator.getSkeletons();
	}

	public ArrayList< HashMap< Integer, Map< String, Object > > > getMeasurementsTimepointList()
	{
		return measurementsTimepointList;
	}

}
