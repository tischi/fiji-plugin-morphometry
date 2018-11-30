package de.embl.cba.morphometry.skeleton;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometrySettings;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;

import java.util.ArrayList;

public class SkeletonAnalyzer< R extends RealType< R > >
{

	final RandomAccessibleInterval< BitType > skeleton;

	public SkeletonAnalyzer( RandomAccessibleInterval< BitType > skeleton )
	{
		this.skeleton = skeleton;
	}


	public SkeletonAnalyzer( ImgLabeling<Integer, IntType> imgLabeling,
							 int label,
							 RandomAccessibleInterval< BitType > skeletons )
	{
		final LabelRegions< Integer > labelRegions = new LabelRegions<>( imgLabeling );

		final LabelRegion< Integer > labelRegion = labelRegions.getLabelRegion( label );

//		Intervals.minAsLongArray( skeletons )


		this.skeleton = null;
	}


	public void run()
	{

//			final RandomAccessibleInterval< BitType > skeletons =
//					Algorithms.createObjectSkeletons(
//							imgLabeling,
//							3, // TODO: Make a parameter
//							settings.opService
//					);

//			this.skeletons.add( skeletons )

	}

	public long getNumBranchPoints()
	{
		return 0;
	}

}
