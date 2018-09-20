#include <QApplication>
#include <QDesktopWidget>

// modification in ParaView\Qt\ApplicationComponents\pqSaveStateReaction.cxx

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void pqSaveStateReaction::saveState()
{

  //...
 
  if (fileDialog.exec() == QDialog::Accepted)
    {
    QString selectedFile = fileDialog.getSelectedFiles()[0];
    if(selectedFile.endsWith(".py"))
      {
      pqSaveStateReaction::savePythonState(selectedFile);
      }
    else
      {
      pqSaveStateReaction::saveState(selectedFile);
      }

    // ### PATCH_BEGIN ###
      
    // wait for save dialog to close properly..
    QTime t;
    t.start();
    while(t.elapsed() < 1000) {
        QApplication::processEvents();
    }

    // grab whole desktop and save..
    QPixmap screenshot;
    WId root = QApplication::desktop()->winId();
    screenshot = QPixmap::grabWindow(root);

    QWidget *w = QApplication::activeWindow();
    if(w) {
        QPoint pos = w->mapToGlobal(QPoint(0,0));
        QPixmap buf = screenshot.copy(pos.x(), pos.y(), w->width(), w->height());
        screenshot = buf;
    }

    QString screenshotName = selectedFile;
    int index = selectedFile.lastIndexOf('.');
    if( index > 0)
        screenshotName = selectedFile.left( index);
        screenshot.save(QString(screenshotName+".png"));
    }
    
    // ### PATCH_END ###

}
