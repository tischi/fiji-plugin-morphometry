package de.embl.cba.morphometry.splitting;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.ObjectMeasurements;
import de.embl.cba.morphometry.microglia.MicrogliaSettings;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegionCursor;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;

import java.util.ArrayList;
import java.util.HashMap;

import static de.embl.cba.morphometry.tracking.TrackingUtils.computeOverlaps;

public class TrackingSplitter< T extends RealType< T > & NativeType< T > >
{

	final ArrayList< RandomAccessibleInterval< BitType > > masks;
	final ArrayList< RandomAccessibleInterval< T > > intensities;

	private int nextId;

	private ArrayList< RandomAccessibleInterval< BitType > > splitMasks;
	final MicrogliaSettings settings;


	public TrackingSplitter( ArrayList< RandomAccessibleInterval< BitType > > masks,
							 ArrayList< RandomAccessibleInterval< T > > intensities,
							 MicrogliaSettings settings )
	{
		this.masks = masks;
		this.intensities = intensities;
		this.settings = settings;
	}


	public void run()
	{

		int tMin = 0;  // at this point the movie is already cropped in time, such that we can process the full movie
		int tMax = masks.size() - 1;
		int t = tMin;

		splitMasks = new ArrayList<>( );


		/**
		 * Process first time-point
		 */

		Utils.log( "\nRunning ShapeAndIntensitySplitter on frame " + settings.tMin );
		final ShapeAndIntensitySplitter splitter = new ShapeAndIntensitySplitter( masks.get( t ), intensities.get( t ), settings );
		splitter.run();
		splitMasks.add( splitter.getSplitMask() );

		nextId = Utils.getNumObjects( splitMasks.get( tMin ) );

		boolean showSplittingAttempts = false;

		/**
		 * Process subsequent time-points
		 */


		RandomAccessibleInterval< IntType > previousLabeling = Utils.asImgLabeling( splitMasks.get( tMin ) ).getSource();

		for ( t = tMin + 1; t <= tMax; ++t )
		{

			Utils.log( "\nProcessing frame " + ( t + 1 ) );
			final ImgLabeling< Integer, IntType > currentImgLabeling = Utils.asImgLabeling( masks.get( t ) );
			RandomAccessibleInterval< IntType > currentLabeling = currentImgLabeling.getSource();

			HashMap< Integer, Integer > numObjectsPerRegion = new HashMap<>(  );

			LabelRegions< Integer > labelRegions = new LabelRegions( currentImgLabeling );

			for ( LabelRegion< Integer > region : labelRegions )
			{
				final HashMap< Integer, Long > overlaps = computeOverlaps( previousLabeling, region );

				int numObjectsInRegion = 1;

				if ( overlaps.size() > 1 )
				{

					Utils.log( "Object at "
							+ " x = " + (int) region.getCenterOfMass().getDoublePosition( 0 )
							+ " y = " + (int) region.getCenterOfMass().getDoublePosition( 1 )
							+ " overlaps with " + overlaps.size() + " objects in previous frame." );

					if ( overlaps.size() == 2 )
					{
						boolean splitObjects = isReallyTwoObjects( t, previousLabeling, currentLabeling, region, overlaps );

						if ( splitObjects )
						{
							numObjectsInRegion = 2;
						}
					}
				}

				numObjectsPerRegion.put( region.getLabel(), numObjectsInRegion );

			}

			final RandomAccessibleInterval< BitType > splitMask = Utils.copyAsArrayImg( masks.get( t ) );

			Algorithms.splitTouchingObjects(
					currentImgLabeling,
					intensities.get( t ),
					splitMask,
					numObjectsPerRegion,
					( int ) ( settings.minimalObjectCenterDistance / settings.workingVoxelSize ),
					( long ) ( settings.minimalObjectSize / Math.pow( settings.workingVoxelSize, splitMask.numDimensions() ) ),
					( int ) ( settings.maximalWatershedLength / settings.workingVoxelSize ),
					settings.opService,
					true,
					false);

			splitMasks.add( splitMask );

			previousLabeling = Utils.asImgLabeling( splitMask ).getSource();
		}
	}

	public boolean isReallyTwoObjects( int t, RandomAccessibleInterval< IntType > previousLabeling, RandomAccessibleInterval< IntType > currentLabeling, LabelRegion< Integer > region, HashMap< Integer, Long > overlaps )
	{

		boolean splitObjects =  true;

		final double currentIntensity = ObjectMeasurements.measureBgCorrectedSumIntensity(
				currentLabeling,
				region.getLabel(),
				intensities.get( t ) );

		Utils.log( "Object intensity: " + (long) currentIntensity );

		final HashMap< Integer, Double > previousIntensities = new HashMap<>();

		for ( int label : overlaps.keySet() )
		{
			final double intensity = ObjectMeasurements.measureBgCorrectedSumIntensity(
					previousLabeling,
					label,
					intensities.get( t ) );

			previousIntensities.put( label , intensity );

			final double overlapFraction = 1.0 * overlaps.get( label ).longValue() / region.size();
			Utils.log( "Previous object intensity: " + (long) intensity );
			Utils.log( "Overlap pixel fraction: " + overlapFraction );

			if ( overlapFraction < settings.minimalOverlapFraction ) splitObjects = false;

		}

		double previousIntensitySum = getPreviousIntensitySum( previousIntensities );

		final double sumIntensityRatio = currentIntensity / previousIntensitySum;

		Utils.log( "Intensity ratio: " + sumIntensityRatio );

		if ( sumIntensityRatio < settings.minimalSumIntensityRatio ) splitObjects = false;
		if ( sumIntensityRatio > settings.maximalSumIntensityRatio ) splitObjects = false;

		Utils.log( "Split objects: " + splitObjects );

		return splitObjects;
	}

	private double getPreviousIntensitySum( HashMap< Integer, Double > previousIntensities )
	{
		double previousIntensitySum = 0;

		for ( double intensity : previousIntensities.values() )
		{
			previousIntensitySum += intensity;

		}
		return previousIntensitySum;
	}

	public ArrayList< RandomAccessibleInterval< BitType > > getSplitMasks()
	{
		return splitMasks;
	}



	public ArrayList< RandomAccessibleInterval< IntType > > initUpdatedLabelings( int tMin )
	{
//		final ArrayList< RandomAccessibleInterval< IntType > > updatedLabelings = new ArrayList<>();
//		updatedLabelings.add( masks.get( tMin ).getSource() );
//		return updatedLabelings;
		return null;
	}

	public int computeObjectId( HashMap< Integer, Long > overlaps )
	{
		int objectId;

		if( overlaps.size() == 0 )
		{
			objectId = nextId++;
		}
		else
		{
			objectId = getMaxOverlapLabel( overlaps );
		}

		return objectId;
	}

	public int getMaxOverlapLabel( HashMap< Integer, Long > overlaps )
	{
		int maxOverlapLabel = 0;
		long maxOverlapCount = Long.MIN_VALUE;

		for ( Integer label : overlaps.keySet() )
		{
			final long overlapCount = overlaps.get( label ).longValue();

			if ( overlapCount > maxOverlapCount )
			{
				maxOverlapCount = overlapCount;
				maxOverlapLabel = label;
			}
		}
		return maxOverlapLabel;
	}

	public HashMap< Integer, Long > computeOverlaps(
			RandomAccessibleInterval< IntType > previousLabeling,
			LabelRegion region )
	{
		final HashMap< Integer, Long > overlaps = new HashMap<>();

		final RandomAccess< IntType > previousLabelsAccess = previousLabeling.randomAccess();
		final LabelRegionCursor cursor = region.cursor();

		while ( cursor.hasNext() )
		{
			cursor.fwd();
			previousLabelsAccess.setPosition( cursor );

			final int label = previousLabelsAccess.get().getInteger();

			if ( label != 0 )
			{
				addOverlap( overlaps, label );
			}
		}

		return overlaps;
	}

	public void addOverlap( HashMap< Integer, Long > overlaps, int integer )
	{
		if ( overlaps.keySet().contains( integer ) )
		{
			overlaps.put( integer,  overlaps.get( integer ).longValue() + 1 );
		}
		else
		{
			overlaps.put( integer,  1L );
		}
	}


}
