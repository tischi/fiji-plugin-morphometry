package de.embl.cba.morphometry.skeleton;

import net.imagej.ops.OpService;
import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

public class SkeletonAnalyzer< R extends RealType< R > >
{

	final RandomAccessibleInterval< BitType > skeleton;
	final OpService opService;

	double skeletonLength;
	long numBranchPoints;
	private RandomAccessibleInterval<BitType> branchpoints;

	public SkeletonAnalyzer( RandomAccessibleInterval< BitType > skeleton, OpService opService )
	{
		this.skeleton = skeleton;
		this.opService = opService;

		run();
	}


	private void run()
	{
		measureLength();

		branchpoints = measureBranchpoints();

	}

	private RandomAccessibleInterval< BitType > measureBranchpoints()
	{
		RandomAccessibleInterval< BitType > branchpoints = Skeletons.branchpoints( skeleton );

		numBranchPoints = measureSum( branchpoints );

		return branchpoints;
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
