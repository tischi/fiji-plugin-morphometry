package de.embl.cba.morphometry.measurements;

import bdv.util.RandomAccessibleIntervalSource;
import de.embl.cba.bdv.utils.argbconversion.SelectableRealVolatileARGBConverter;
import de.embl.cba.bdv.utils.argbconversion.VolatileARGBConvertedRealSource;
import de.embl.cba.tables.TableUtils;
import ij.IJ;
import ij.ImagePlus;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.util.Util;
import net.imglib2.view.Views;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;

import javax.swing.*;
import java.io.File;
import java.io.IOException;

public class ObjectMeasurementsReviewCommand implements Command
{

	@Parameter ( label = "Results table" )
	public File inputTableFile;

	@Parameter ( label = "Label masks (single channel, 2D or 3D)" )
	public File inputLabelMasksFile;

	@Parameter ( label = "Intensities (optional)")
	public File inputIntensitiesFile;

	@Override
	public void run()
	{
		final JTable jTable = loadTable( inputTableFile );

		final RandomAccessibleIntervalSource labels = loadImagesAsRAISource( inputLabelMasksFile );
		final SelectableRealVolatileARGBConverter argbConverter = new SelectableRealVolatileARGBConverter();
		final VolatileARGBConvertedRealSource argbSource = new VolatileARGBConvertedRealSource( raiSource,  argbConverter );

	}

	public JTable loadTable( File file )
	{
		try
		{
			return TableUtils.loadTable( file, "\t" );
		}
		catch ( IOException e )
		{
			e.printStackTrace();
		}

		return null;
	}

	public static RandomAccessibleIntervalSource loadImagesAsRAISource( File file )
	{
		final ImagePlus imagePlus = IJ.openImage( file.toString() );

		RandomAccessibleInterval< RealType > wrap = ImageJFunctions.wrapReal( imagePlus );

		if ( imagePlus.getNSlices() == 1 )
		{
			// needs to be at least 3D
			wrap = Views.addDimension( wrap, 0, 0 );
		}

		return new RandomAccessibleIntervalSource( wrap, Util.getTypeFromInterval( wrap ), imagePlus.getTitle() );
	}
}
