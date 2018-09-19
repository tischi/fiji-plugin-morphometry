package de.embl.cba.morphometry.tracking;

import net.imglib2.RandomAccessibleInterval;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.numeric.integer.IntType;

public class Tracker
{

	final RandomAccessibleInterval< IntType > labeling;

	public Tracker( RandomAccessibleInterval< IntType > labeling )
	{
		this.labeling = labeling;
	}


	public void run()
	{
		int tMax = 0;

	}
}
