package de.embl.cba.morphometry.skeleton;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometrySettings;
import de.embl.cba.morphometry.microglia.MicrogliaTrackingSettings;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.view.Views;

import java.util.ArrayList;

public class SkeletonCreator< T extends RealType< T > & NativeType< T > >
{

	final ArrayList< RandomAccessibleInterval< BitType > > masks;
	final MicrogliaMorphometrySettings settings;

	private ArrayList< RandomAccessibleInterval< BitType > > skeletons;

	public SkeletonCreator( ArrayList< RandomAccessibleInterval< BitType > > masks,
							MicrogliaMorphometrySettings settings )
	{
		this.masks = masks;
		this.settings = settings;
	}

	public void run()
	{

		int tMin = 0;  // at this point the movie is already cropped in time, such that we can process the full movie
		int tMax = masks.size() - 1;

		skeletons = new ArrayList<>( );

		for ( int t = tMin; t <= tMax; ++t )
		{
			Utils.log( "Computing skeletons for frame " + ( t + 1 ) );

			final ImgLabeling< Integer, IntType > imgLabeling = Utils.asImgLabeling( masks.get( t ) );

			final RandomAccessibleInterval< BitType > skeletons =
					Algorithms.createObjectSkeletons(
						imgLabeling,
						3, // TODO: Make a parameter
						settings.opService
					);

			this.skeletons.add( skeletons );
		}

	}

	public ArrayList< RandomAccessibleInterval< BitType > > getSkeletons()
	{
		return skeletons;
	}

}
