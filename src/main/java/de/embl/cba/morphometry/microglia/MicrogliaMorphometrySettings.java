package de.embl.cba.morphometry.microglia;

import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.io.File;

public class MicrogliaMorphometrySettings<T extends RealType<T> & NativeType< T > >
{

	public OpService opService;
	public long labelMapChannelIndex = 2;

	public double[] inputCalibration;
	public boolean showIntermediateResults = false;
}
