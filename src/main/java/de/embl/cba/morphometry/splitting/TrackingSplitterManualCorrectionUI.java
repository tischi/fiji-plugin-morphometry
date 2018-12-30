package de.embl.cba.morphometry.splitting;

import de.embl.cba.morphometry.SyncWindowsHack;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.regions.Regions;
import de.embl.cba.morphometry.tracking.MaximalOverlapTracker;
import ij.IJ;
import ij.ImagePlus;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

public class TrackingSplitterManualCorrectionUI < T extends RealType< T > & NativeType< T > > extends JPanel
{
	private final long minimalObjectSizeInPixels;
	private JFrame frame;
	private boolean isFinished;
	private ImagePlus editedLabelsImp;
	private ArrayList< RandomAccessibleInterval< T > > labels;
	private SyncWindowsHack syncWindows;
	private static Point frameLocation;
	private static Point editedLabelsImpLocation;

	public TrackingSplitterManualCorrectionUI(
			ArrayList< RandomAccessibleInterval< T > > labels,
			long minimalObjectSizeInPixels )
	{
		this.isFinished = false;
		this.minimalObjectSizeInPixels = minimalObjectSizeInPixels;

		showNewLabelsImagePlusForEditing( labels );

		add( updateLabelsButton() );

		add( nextFrameButton() );

		showPanel();
	}

	public void showNewLabelsImagePlusForEditing( ArrayList< RandomAccessibleInterval< T > > labels )
	{
		editedLabelsImp = Utils.labelingsAsImagePlus( labels );
		editedLabelsImp.show();
		if ( editedLabelsImpLocation != null ) editedLabelsImp.getWindow().setLocation( editedLabelsImpLocation );
		editedLabelsImp.setT( editedLabelsImp.getNFrames() );
		editedLabelsImp.updateImage();
		editedLabelsImp.setActivated();
		IJ.run( editedLabelsImp, "Enhance Contrast", "saturated=0.00");
		IJ.run("Brightness/Contrast...");

		// in case other images, e.g. the raw intensities are shown at the same time
		syncWindows = new SyncWindowsHack();
		syncWindows.syncAll();
	}

	public JButton updateLabelsButton()
	{
		final JButton button = new JButton( "Update Labels" );
		button.addActionListener( new ActionListener()
		{
			@Override
			public void actionPerformed( ActionEvent e )
			{
				labels = runMaximalOverlapTrackerOnEditedImagePlus();

				closeCurrentEditedLabelsImagePlus();

				showNewLabelsImagePlusForEditing( labels );
			}
		} );

		return button;
	}

	public void closeCurrentEditedLabelsImagePlus()
	{
		editedLabelsImp.changes = false;
		editedLabelsImpLocation = editedLabelsImp.getWindow().getLocation();
		editedLabelsImp.close();

		syncWindows.close();
	}

	public JButton nextFrameButton()
	{
		final JButton button = new JButton( "Next frame" );
		button.addActionListener( new ActionListener()
		{
			@Override
			public void actionPerformed( ActionEvent e )
			{
				labels = runMaximalOverlapTrackerOnEditedImagePlus();
				isFinished = true;
				closeCurrentEditedLabelsImagePlus();
				frameLocation = frame.getLocation();
				frame.dispose();
			}
		} );
		return button;
	}


	private void showPanel() {

		//Create and set up the window.
		frame = new JFrame("Manual editing");

		//Create and set up the content pane.
		this.setOpaque(true); //content panes must be opaque
		frame.setContentPane(this);

		//Display the window.
		frame.pack();
		if ( frameLocation != null ) frame.setLocation( frameLocation );
		frame.setVisible( true );

	}

	public boolean isFinished()
	{
		return isFinished;
	}

	public ArrayList< RandomAccessibleInterval< T > > runMaximalOverlapTrackerOnEditedImagePlus()
	{
		final ArrayList< RandomAccessibleInterval< T > > labels = Utils.get2DImagePlusMovieAsFrameList( editedLabelsImp, 1 );

		// Due to the editing small unconnected regions of pixels may occur
		Regions.removeSmallRegionsInMasks( labels, minimalObjectSizeInPixels );

		// TODO: if below turns out to be too slow, only do it for the last two frames
		final MaximalOverlapTracker maximalOverlapTracker = new MaximalOverlapTracker( labels );
		maximalOverlapTracker.run();

		return maximalOverlapTracker.getLabelings();
	}

	public ArrayList< RandomAccessibleInterval< T > > getLabels()
	{
		return labels;
	}


}
