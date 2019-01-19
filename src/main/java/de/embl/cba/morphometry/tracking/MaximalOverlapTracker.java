package de.embl.cba.morphometry.tracking;

import de.embl.cba.morphometry.Utils;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.util.Intervals;

import java.util.ArrayList;
import java.util.HashMap;

public class MaximalOverlapTracker< T extends RealType< T > & NativeType< T > >
{

	ArrayList< RandomAccessibleInterval< T > > masks;

	private int nextId;

	private ArrayList< RandomAccessibleInterval< IntType > > labelings;

	public MaximalOverlapTracker( ArrayList< RandomAccessibleInterval< T > > masks )
	{
		this.masks = masks;
	}

	public void run()
	{

		int tMin = 0;
		int tMax = masks.size() - 1;

		int t = tMin;

		nextId = Utils.getNumObjects( masks.get( t ) );

		RandomAccessibleInterval< IntType > previousLabeling = Utils.asImgLabeling( masks.get( tMin ), ConnectedComponents.StructuringElement.FOUR_CONNECTED ).getSource();
		LabelRegions< Integer > labelRegions;

		labelings = initLabelings( masks.get( tMin )  );

		for ( t = tMin + 1; t <= tMax; ++t )
		{
			RandomAccessibleInterval< IntType > currentLabeling = Utils.asImgLabeling( masks.get( t ), ConnectedComponents.StructuringElement.FOUR_CONNECTED ).getSource();

			RandomAccessibleInterval< IntType > updatedLabeling = ArrayImgs.ints( Intervals.dimensionsAsLongArray( currentLabeling ) );

			labelRegions = new LabelRegions( Utils.asImgLabeling( currentLabeling, ConnectedComponents.StructuringElement.FOUR_CONNECTED ) );

			for ( LabelRegion< Integer > region : labelRegions )
			{
				final HashMap< Integer, Long > overlaps = TrackingUtils.computeOverlaps( previousLabeling, region );

				int objectId = TrackingUtils.computeObjectId( overlaps, nextId );

				Utils.drawObject( updatedLabeling, region, objectId );
			}

			labelings.add( updatedLabeling);

			previousLabeling = updatedLabeling;
		}
	}

	public ArrayList< RandomAccessibleInterval< IntType > > getLabelings()
	{
		return labelings;
	}


	public ArrayList< RandomAccessibleInterval< IntType > > initLabelings( RandomAccessibleInterval< T > mask )
	{
		final ArrayList< RandomAccessibleInterval< IntType > > updatedLabelings = new ArrayList<>();
		updatedLabelings.add( Utils.asImgLabeling( mask, ConnectedComponents.StructuringElement.FOUR_CONNECTED ).getIndexImg() );
		return updatedLabelings;
	}

}
