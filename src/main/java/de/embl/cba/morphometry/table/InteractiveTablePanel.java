package de.embl.cba.morphometry.table;

import net.imagej.table.GenericTable;

import javax.swing.*;
import java.awt.*;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;


public class InteractiveTablePanel extends JPanel implements MouseListener, KeyListener {
    private boolean DEBUG = false;
    JTable table;
    JFrame frame;
    JScrollPane scrollPane;

    public InteractiveTablePanel( GenericTable genericTable ) {

        super( new GridLayout(1, 0 ) );

        table = TableUtils.asJTable( genericTable );
        table.setPreferredScrollableViewportSize( new Dimension(500, 200) );
        table.setFillsViewportHeight( true );
        table.setAutoCreateRowSorter( true );
        table.setRowSelectionAllowed( true );
        table.addMouseListener(this );
        table.addKeyListener(this );

        scrollPane = new JScrollPane( table );
        add( scrollPane );

    }

    public void showTable() {

        //Create and set up the window.
        frame = new JFrame("Table");

        //Create and set up the content pane.
        this.setOpaque(true); //content panes must be opaque
        frame.setContentPane(this);

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }

    public void reactToAction() {


        int rs = table.getSelectedRow();
        int r = table.convertRowIndexToModel(rs);
        

//        ImagePlus imp = IJ.getImage();

//        float x = new Float(table.getModel().getValueAt(r, 1).toString());
//        float y = new Float(table.getModel().getValueAt(r, 2).toString());
//        float z = new Float(table.getModel().getValueAt(r, 3).toString());
//        int t = new Integer(table.getModel().getValueAt(r, 4).toString());
//        int id = new Integer(table.getModel().getValueAt(r, 5).toString());
//
//        // Values in the table are one-based already
//        imp.setPosition(0,(int)z, t);
//        Roi pr = new PointRoi(x,y);
//        pr.setPosition(0,(int)z, t);
//        imp.setRoi(pr);

    }

    @Override
    public void mouseClicked(MouseEvent e)
    {
        reactToAction();
    }

    @Override
    public void mousePressed(MouseEvent e) {

    }

    @Override
    public void mouseReleased(MouseEvent e) {

    }

    @Override
    public void mouseEntered(MouseEvent e) {

    }

    @Override
    public void mouseExited(MouseEvent e) {

    }

    @Override
    public void keyTyped(KeyEvent e) {

    }

    @Override
    public void keyPressed(KeyEvent e) {

    }

    @Override
    public void keyReleased(KeyEvent e)
    {
        reactToAction();
    }
}
