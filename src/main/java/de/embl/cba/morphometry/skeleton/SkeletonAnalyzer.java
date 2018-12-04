package de.embl.cba.morphometry.skeleton;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometrySettings;
import net.imagej.ops.OpService;
import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.morphology.table2d.Branchpoints;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegionCursor;
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
	final OpService opService;

	double skeletonLength;
	long numBranchPoints;
	private RandomAccessibleInterval< BitType > branchpoints;

	public SkeletonAnalyzer( RandomAccessibleInterval< BitType > skeleton, OpService opService )
	{
		this.skeleton = skeleton;
		this.opService = opService;

		run();
	}


	private void run()
	{
		measureLength();

		measureBranchpoints();

	}

	private void measureBranchpoints()
	{
		branchpoints = ArrayImgs.bits( Intervals.dimensionsAsLongArray( skeleton ) );

		Branchpoints.branchpoints(
				Views.extendBorder( Views.zeroMin( skeleton ) ),
				Views.flatIterable( branchpoints ) );

		numBranchPoints = measureSum( branchpoints );
	}

	public void measureLength()
	{
		skeletonLength = measureSum( skeleton );
	}

	public long getNumBranchPoints()
	{
		return numBranchPoints;
	}

	public RandomAccessibleInterval< BitType > getBranchpoints()
	{
		return branchpoints;
	}


	public double getSkeletonLength() { return skeletonLength; }

	public static long measureSum( RandomAccessibleInterval< BitType > rai )
	{
		final Cursor< BitType > cursor = Views.iterable( rai ).cursor();

		long sum = 0;

		while ( cursor.hasNext() )
		{
			sum += cursor.next().getRealDouble();
		}

		return sum;
	}



}
