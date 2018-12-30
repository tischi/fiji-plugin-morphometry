package de.embl.cba.morphometry.tracking;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.roi.labeling.*;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.util.Intervals;

import java.util.ArrayList;
import java.util.HashMap;

public class Tracker < T extends RealType< T > & NativeType< T > >
{

	final ArrayList< ImgLabeling< Integer, IntType > > imgLabelings;
	final ArrayList< RandomAccessibleInterval< T > > intensities;

	private int nextId;

	private ArrayList< RandomAccessibleInterval< IntType > > updatedLabelings;


	public Tracker( ArrayList< ImgLabeling< Integer, IntType > > imgLabelings,
					ArrayList< RandomAccessibleInterval< T > > intensities )
	{
		this.imgLabelings = imgLabelings;
		this.intensities = intensities;
	}


	public void run()
	{

		int tMin = 0;
		int tMax = imgLabelings.size() - 1;

		int t = tMin;

		nextId = getNumObjects( t );

		RandomAccessibleInterval< IntType > previousLabeling = imgLabelings.get( tMin ).getSource();
		LabelRegions< Integer > labelRegions;

		updatedLabelings = initUpdatedLabelings( tMin );

		for ( t = tMin + 1; t <= tMax; ++t )
		{
			RandomAccessibleInterval< IntType > currentLabeling = imgLabelings.get( t ).getSource();;
			RandomAccessibleInterval< IntType > updatedLabeling = ArrayImgs.ints( Intervals.dimensionsAsLongArray( currentLabeling ) );

			labelRegions = new LabelRegions( Utils.asImgLabeling( currentLabeling ) );

			for ( LabelRegion< Integer > region : labelRegions )
			{
				final HashMap< Integer, Long > overlaps = computeOverlaps( previousLabeling, region );

				if ( overlaps.size() == 2 )
				{

					Logger.log( "Overlap with two objects" );

					Logger.log( "Time point (one based): " + ( t + 1 ) );
					Logger.log( "Object label: " + region.getLabel() );

					final double currentIntensity = Measurements.measureBgCorrectedSumIntensity(
							currentLabeling,
							region.getLabel(),
							intensities.get( t ) );

					Logger.log( "Object intensity: " + (long) currentIntensity );

					final HashMap< Integer, Double > previousIntensities = new HashMap<>();

					for ( int label : overlaps.keySet() )
					{
						final double intensity = Measurements.measureBgCorrectedSumIntensity(
								previousLabeling,
								label,
								intensities.get( t ) );

						previousIntensities.put( label , intensity );

						Logger.log( "Previous object intensity: " + (long) intensity );
						Logger.log( "Overlap pixel fraction: " + 1.0 * overlaps.get( label ).longValue() / region.size() );
					}

					double previousIntensitySum = getPreviousIntensitySum( previousIntensities );

					Logger.log( "Intensity ratio: " +  currentIntensity / previousIntensitySum );



				}

				int objectId = computeObjectId( overlaps );

				Utils.drawObject( updatedLabeling, region, objectId );
			}

			updatedLabelings.add( updatedLabeling);

			previousLabeling = updatedLabeling;
		}
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

	public ArrayList< RandomAccessibleInterval< IntType > > getUpdatedLabelings()
	{
		return updatedLabelings;
	}



	public ArrayList< RandomAccessibleInterval< IntType > > initUpdatedLabelings( int tMin )
	{
		final ArrayList< RandomAccessibleInterval< IntType > > updatedLabelings = new ArrayList<>();
		updatedLabelings.add( imgLabelings.get( tMin ).getSource() );
		return updatedLabelings;
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
		final LabelRegions labelRegions = new LabelRegions( imgLabelings.get( t ) );
		return labelRegions.getExistingLabels().size() - 1;
	}
}
