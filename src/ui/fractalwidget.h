#ifndef FRACTALWIDGET_H
#define FRACTALWIDGET_H

#include "config.h"
#include <QWidget>
#ifdef USE_OPENGL
#include <QOpenGLWidget>
#endif
class QImage;
class QPoint;

#ifdef USE_OPENGL
class FractalWidget : public QOpenGLWidget
#else
class FractalWidget : public QWidget
#endif
{
    Q_OBJECT
  private:
    struct image *m_image = NULL;
    QSize m_sizeHint;
    QPointF m_mousePosition = QPointF(0.0, 0.0);
    /* Selection zoom. The rectangle lives here, in device pixels like
     * everything the engine is told about, and is drawn by this widget rather
     * than into the fractal image: an overlay in the image would have to be
     * saved and restored around every frame, which is the machinery that
     * smears text across the picture when it goes wrong. */
    struct uih_context *m_uih = NULL;
    QPointF m_selectionStart;
    bool m_selecting = false;

  protected:
    void mouseMoveEvent(QMouseEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);
#ifdef USE_OPENGL
    void paintGL();
    void resizeGL(int w, int h);
#else
    void paintEvent(QPaintEvent *event);
#endif
  public:
    FractalWidget();
    QSize sizeHint() const;
    QPointF mousePosition();
    // The size the fractal image has to be calculated at, in real device
    // pixels rather than the logical ones the widget is laid out in. The two
    // differ by the display scaling factor, which is what Windows changes when
    // the text size is raised above 100%.
    QSize imageSize() const;
    void setImage(struct image *image);
    void setContext(struct uih_context *uih) { m_uih = uih; }
    /* The dragged rectangle, forced to the widget's aspect ratio, in device
     * pixels; null while nothing is being dragged. */
    QRectF selection() const;
};

#endif // FRACTALWIDGET_H
