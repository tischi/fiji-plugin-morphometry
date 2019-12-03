package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.fccf.FCCF;
import de.embl.cba.tables.Tables;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import loci.common.DebugTools;
import org.scijava.command.Command;
import org.scijava.command.CommandService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.widget.Button;

import javax.swing.*;
import javax.swing.table.TableColumn;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

@Plugin(type = Command.class, menuPath = "Plugins>EMBL>FCCF>BD View Images From Table" )
public class BDOpenTableCommand implements Command
{
	@Parameter
	public LogService logService;

	@Parameter
	public CommandService commandService;

	@Parameter ( label = "Image Table" )
	public File imageTablePath;

	public void run()
	{
		final JTable jTable = loadTable();

		jTable.getColumnModel().getColumnIndex( "class" );
//		Tables.columnMin(  )
//		final TableColumn column = jTable.getColumn();

		final ArrayList< String > classes = new ArrayList<>();
		classes.add( "Hello" );
		classes.add( "World" );
		BDImageViewingCommand.classChoices = classes;
		commandService.run( BDImageViewingCommand.class, true );
	}

	public JTable loadTable()
	{
		// TODO: load table
		final long currentTimeMillis = System.currentTimeMillis();
		IJ.log("Loading table; please wait...");
		final JTable jTable = Tables.loadTable( " " );
		IJ.log( "Loaded table in " + ( System.currentTimeMillis() - currentTimeMillis ) + " ms." );
		return jTable;
	}


}
