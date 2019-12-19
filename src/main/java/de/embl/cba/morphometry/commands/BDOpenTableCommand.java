package de.embl.cba.morphometry.commands;

import de.embl.cba.tables.Tables;
import ij.IJ;
import org.scijava.command.Command;
import org.scijava.command.CommandService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.widget.Button;

import javax.swing.*;
import java.io.File;
import java.util.List;

@Plugin(type = Command.class, menuPath = "Plugins>EMBL>FCCF>BD View Images from Table" )
public class BDOpenTableCommand implements Command
{
	@Parameter
	public LogService logService;

	@Parameter
	public CommandService commandService;

	@Parameter ( label = "Image Table" )
	public File imageTableFile;

	@Parameter ( label = "Glimpse Table", callback = "logTableInfo" )
	public Button logTableInfo;

	@Parameter ( label = "Image Path Column Name" )
	public String imagePathColumnName = "path";

	@Parameter ( label = "Gate Column Name" )
	public String gateColumnName = "Gate";

	private JTable jTable;

	private String recentImageTablePath = "";

	public void run()
	{
		loadTable();

		BDImageViewingCommand.jTable = jTable;
		BDImageViewingCommand.imagesRootDir = imageTableFile.getParent();
		BDImageViewingCommand.imagePathColumnName = imagePathColumnName;
		BDImageViewingCommand.gateColumnName = gateColumnName;

		commandService.run( BDImageViewingCommand.class, true );
	}

	private void loadTable()
	{
		if ( ! imageTableFile.exists() )
		{
			logService.error( "Table file does not exist: " + imageTableFile );
			throw new UnsupportedOperationException( "Could not open file: " + imageTableFile);
		}

		if ( recentImageTablePath.equals( imageTableFile.getAbsolutePath() ) ) return;

		final long currentTimeMillis = System.currentTimeMillis();
		IJ.log("Loading table; please wait...");
		jTable = Tables.loadTable( imageTableFile.getAbsolutePath() );
		IJ.log( "Loaded table in " + ( System.currentTimeMillis() - currentTimeMillis ) + " ms." );

		recentImageTablePath = imageTableFile.getAbsolutePath();
	}

	public void logTableInfo()
	{
		loadTable();

		IJ.log( "# Table Info"  );
		IJ.log( "Number of rows: " + jTable.getRowCount() );
		final List< String > columnNames = Tables.getColumnNames( jTable );
		for ( String columnName : columnNames )
		{
			final int columnIndex = jTable.getColumnModel().getColumnIndex( columnName );

			String firstRows = "";
			for ( int i = 0; i < 5; i++ )
			{
				firstRows += jTable.getValueAt( 0, columnIndex );
				firstRows += ", ";
			}
			firstRows += "...";

			IJ.log( columnName + ": " + firstRows );
		}
	}


}
