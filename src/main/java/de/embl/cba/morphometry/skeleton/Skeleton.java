package de.embl.cba.morphometry.skeleton;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.ObjectMeasurements;
import de.embl.cba.morphometry.microglia.MicrogliaSettings;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.NonBlockingGenericDialog;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
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

public class Skeleton< T extends RealType< T > & NativeType< T > >
{

	final ArrayList< RandomAccessibleInterval< BitType > > masks;
	final MicrogliaSettings settings;

	private ArrayList< RandomAccessibleInterval< BitType > > skeletons;

	public Skeleton( ArrayList< RandomAccessibleInterval< BitType > > masks,
					 MicrogliaSettings settings )
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

			Utils.log( "\nComputing skeletons for frame " + ( t + 1 ) );

			final ImgLabeling< Integer, IntType > imgLabeling = Utils.asImgLabeling( masks.get( t ) );

			final RandomAccessibleInterval< BitType > skeletons =
					Algorithms.createObjectSkeletons(
						imgLabeling,
						3,
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
