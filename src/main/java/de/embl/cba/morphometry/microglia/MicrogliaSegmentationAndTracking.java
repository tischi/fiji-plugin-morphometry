package de.embl.cba.morphometry.microglia;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.segmentation.MicrogliaSegmenter;
import de.embl.cba.morphometry.tracking.SemiAutomatedTrackingSplitter;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.util.ArrayList;

public class MicrogliaSegmentationAndTracking< T extends RealType< T > & NativeType< T > >
{
	private final MicrogliaSettings settings;
	private final ArrayList< RandomAccessibleInterval< T > > intensities;
	private ArrayList< RandomAccessibleInterval< T > > labelings;

	public MicrogliaSegmentationAndTracking(
			ArrayList< RandomAccessibleInterval< T  > > intensities,
			MicrogliaSettings settings )
	{
		this.intensities = intensities;
		this.settings = MicrogliaSettings.addMissingSettings( settings );
	}

	public void run()
	{
		ArrayList< RandomAccessibleInterval< T > > masks = createBinaryMasks( intensities );
		labelings = splitTouchingObjectsAndTrack( intensities, masks );
	}

	private ArrayList< RandomAccessibleInterval< T > > createBinaryMasks(
			ArrayList< RandomAccessibleInterval< T > > intensities )
	{
		ArrayList<  RandomAccessibleInterval< T > > masks = new ArrayList<>();
		for ( long t = 0; t < intensities.size() ; ++t )
		{
			Logger.log("Creating mask for frame " + ( t + 1 ) );
			final MicrogliaSegmenter microgliaSegmenter =
					new MicrogliaSegmenter( intensities.get( ( int ) t ), settings );
			microgliaSegmenter.run();
			masks.add( microgliaSegmenter.getMask() );
		}
		return masks;
	}

	private ArrayList< RandomAccessibleInterval< T > > splitTouchingObjectsAndTrack(
			ArrayList< RandomAccessibleInterval< T > > intensities,
			ArrayList< RandomAccessibleInterval< T > > masks )
	{
		final SemiAutomatedTrackingSplitter splitter =
				new SemiAutomatedTrackingSplitter( masks, intensities, settings );

		if ( labelings != null )
		{
			splitter.setLabelings( labelings );
		}

		splitter.run();

		return splitter.getLabelings();
	}

	public ArrayList< RandomAccessibleInterval< T > > getLabelings( )
	{
		return labelings;
	}

	public void setLabelings( ArrayList< RandomAccessibleInterval < T > > labelings )
	{
		this.labelings = labelings;
	}
}
