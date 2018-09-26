package de.embl.cba.morphometry.splitting;

import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometrySettings;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegionCursor;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;

import java.util.ArrayList;
import java.util.HashMap;

public class TrackingSplitter< T extends RealType< T > & NativeType< T > >
{

	final ArrayList< RandomAccessibleInterval< BitType > > masks;
	final ArrayList< RandomAccessibleInterval< T > > intensities;

	private int nextId;

	private ArrayList< RandomAccessibleInterval< BitType > > splitMasks;
	final MicrogliaMorphometrySettings settings;


	public TrackingSplitter( ArrayList< RandomAccessibleInterval< BitType > > masks,
							 ArrayList< RandomAccessibleInterval< T > > intensities,
							 MicrogliaMorphometrySettings settings )
	{
		this.masks = masks;
		this.intensities = intensities;
		this.settings = settings;
	}


	public void run()
	{

		int tMin = 0;
		int tMax = masks.size() - 1;
		int t = tMin;


		/**
		 * Process first time-point
		 */

		final ShapeAndIntensitySplitter splitter = new ShapeAndIntensitySplitter( masks.get( t ), intensities.get( t ), settings );
		splitter.run();
		final RandomAccessibleInterval splitMask = splitter.getSplitMask();

		ImageJFunctions.show( splitMask, "Split Mask - TimePoint 1");


//
//		nextId = getNumObjects( t );
//
//		RandomAccessibleInterval< IntType > previousLabeling = masks.get( tMin ).getSource();
//		LabelRegions< Integer > labelRegions;
//
//		splitMasks = initUpdatedLabelings( tMin );
//
//		for ( t = tMin + 1; t <= tMax; ++t )
//		{
//			RandomAccessibleInterval< IntType > currentLabeling = masks.get( t ).getSource();;
//			RandomAccessibleInterval< IntType > updatedLabeling = ArrayImgs.ints( Intervals.dimensionsAsLongArray( currentLabeling ) );
//
//			labelRegions = new LabelRegions( Utils.asImgLabeling( currentLabeling ) );
//
//			for ( LabelRegion< Integer > region : labelRegions )
//			{
//				final HashMap< Integer, Long > overlaps = computeOverlaps( previousLabeling, region );
//
//				if ( overlaps.size() == 2 )
//				{
//
//					Utils.log( "Overlap with two objects" );
//
//					Utils.log( "Time point (one based): " + ( t + 1 ) );
//					Utils.log( "Object label: " + region.getLabel() );
//
//					final double currentIntensity = ObjectMeasurements.measureBgCorrectedSumIntensity(
//							currentLabeling,
//							region.getLabel(),
//							intensities.get( t ) );
//
//					Utils.log( "Object intensity: " + (long) currentIntensity );
//
//					final HashMap< Integer, Double > previousIntensities = new HashMap<>();
//
//					for ( int label : overlaps.keySet() )
//					{
//						final double intensity = ObjectMeasurements.measureBgCorrectedSumIntensity(
//								previousLabeling,
//								label,
//								intensities.get( t ) );
//
//						previousIntensities.put( label , intensity );
//
//						Utils.log( "Previous object intensity: " + (long) intensity );
//						Utils.log( "Overlap pixel fraction: " + 1.0 * overlaps.get( label ).longValue() / region.size() );
//					}
//
//					double previousIntensitySum = getPreviousIntensitySum( previousIntensities );
//
//					Utils.log( "Intensity ratio: " +  currentIntensity / previousIntensitySum );
//
//
//
//				}
//
//				int objectId = computeObjectId( overlaps );
//
//				Utils.drawObject( updatedLabeling, region, objectId );
//			}
//
//			splitMasks.add( updatedLabeling);
//
//			previousLabeling = updatedLabeling;
//		}
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

	public int getNumObjects( int t )
	{
		final LabelRegions labelRegions = new LabelRegions( masks.get( t ) );
		return labelRegions.getExistingLabels().size() - 1;
	}
}
