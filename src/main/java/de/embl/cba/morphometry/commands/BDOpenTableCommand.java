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

@Plugin(type = Command.class, menuPath = "Plugins>EMBL>FCCF>BD View Images From Table" )
public class BDOpenTableCommand implements Command
{
	@Parameter
	public LogService logService;

	@Parameter
	public CommandService commandService;

	@Parameter ( label = "Image Table" )
	public File imageTablePath;

	@Parameter ( label = "Load Table and Print Column Names", callback = "printColumnNames" )
	public Button printColumnNames;

	@Parameter ( label = "Image Path Column Name" )
	public String imagePathColumnName = "path";

	@Parameter ( label = "Object Class Column Name" )
	public String objectClassColumnName = "Gate";

	private JTable jTable;

	public void run()
	{
		jTable = loadTable();

		BDImageViewingCommand.jTable = jTable;
		BDImageViewingCommand.imagesRootDir = imageTablePath.getParent();
		BDImageViewingCommand.imagePathColumnName = imagePathColumnName;
		BDImageViewingCommand.gateColumnName = objectClassColumnName;

		commandService.run( BDImageViewingCommand.class, true );
	}

	private JTable loadTable()
	{
		if ( jTable != null ) return jTable;

		final long currentTimeMillis = System.currentTimeMillis();
		IJ.log("Loading table; please wait...");
		final JTable jTable = Tables.loadTable( imageTablePath.getAbsolutePath() );
		IJ.log( "Loaded table in " + ( System.currentTimeMillis() - currentTimeMillis ) + " ms." );
		return jTable;
	}

	public void printColumnNames()
	{
		jTable = loadTable();
		final List< String > columnNames = Tables.getColumnNames( jTable );
		for ( String name : columnNames )
		{
			IJ.log(name );
		}
	}


}
