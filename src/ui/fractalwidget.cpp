#include "fractalwidget.h"

#include <QtGui>
#ifdef USE_OPENGL
#include <QtOpenGL>
#endif

#include "ui.h"
#include "filter.h"
#include "ui_helper.h"

FractalWidget::FractalWidget()
{
    m_image = NULL;
    setMouseTracking(true);
    setAutoFillBackground(false);
    setAttribute(Qt::WA_OpaquePaintEvent, true);
}

// Both of these convert from the logical coordinates Qt lays the widget out
// in to the device pixels the fractal engine works in. They are the same thing
// only while the display scaling is 100%; above that the engine has to be told
// the real pixel count, or it computes a smaller image that then gets stretched
// over the window and comes out blocky.
QPointF FractalWidget::mousePosition()
{
    return m_mousePosition * devicePixelRatioF();
}

QSize FractalWidget::imageSize() const
{
    const qreal ratio = devicePixelRatioF();
    return QSize(qMax(1, qRound(width() * ratio)),
                 qMax(1, qRound(height() * ratio)));
}

void FractalWidget::setImage(struct image *image) { m_image = image; }

QSize FractalWidget::sizeHint() const { return QSize(800, 600); }

#ifdef USE_OPENGL
void FractalWidget::paintGL()
{
    if (m_image) {
        QImage *qimage =
            reinterpret_cast<QImage **>(m_image->data)[m_image->currimage];
        // QImage glimage = QGLWidget::convertToOpenGLFormat(*qimage);
        // glDrawPixels(glimage.width(), glimage.height(), GL_RGBA,
        //              GL_UNSIGNED_BYTE, glimage.bits());
        // For some reason, Qt 6 requires mirroring the image and using another color space.
        // The old convertToOpenGLFormat is no longer supported in Qt 6.
        QImage mirrored = qimage->mirrored(false, true);
        glDrawPixels(qimage->width(), qimage->height(), GL_BGRA_EXT,
                  GL_UNSIGNED_BYTE, mirrored.bits());
    }
}

void FractalWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, w, 0, h, -1, 1);
    glMatrixMode(GL_MODELVIEW);
}
#else
void FractalWidget::paintEvent(QPaintEvent */*event*/)
{
    if (m_image) {
        QPainter painter(this);
        QImage *qimage =
            reinterpret_cast<QImage **>(m_image->data)[m_image->currimage];
        painter.setCompositionMode(QPainter::CompositionMode_Source);
        // The image is calculated in device pixels (see imageSize), so drawing
        // it across the widget's logical rectangle lands it one image pixel
        // per device pixel, with no resampling.
        //
        // Not by setting a device pixel ratio on the image, which is how this
        // was first done: that is not a passive label. It changes the unit
        // every later QPainter on the image works in, and the overlays -- the
        // status bar, the messages -- are drawn through exactly such a
        // painter, in xprint(). They ended up placed and sized in logical
        // units while the engine and the save/restore in wstack.cpp address
        // the buffer by raw pixel, so at 125% scaling the text was drawn a
        // quarter further right and a quarter wider than the area saved for
        // it. Whatever fell outside was never restored, and the zoom smeared
        // it across the image.
        painter.drawImage(QRectF(0, 0, width(), height()), *qimage);

        /* The selection goes on top of the image and nowhere near it: drawn
         * here every paint, it leaves no trace to be cleaned up and cannot
         * reach the buffers the zoom engine reuses. */
        QRectF r = selection();
        if (!r.isNull()) {
            painter.setCompositionMode(QPainter::CompositionMode_SourceOver);
            painter.setPen(QPen(Qt::white, 1, Qt::DashLine));
            painter.drawRect(r);
            painter.setPen(QPen(Qt::black, 1, Qt::DashLine));
            painter.drawRect(r.adjusted(1, 1, -1, -1));
        }
    }
}
#endif

/* Forced to the widget's aspect ratio by growing the short side, never by
 * cropping the long one: a selection is a promise about what will be on
 * screen, and trimming it would quietly break that promise. The engine is told
 * the same rectangle, so what is drawn here and what is zoomed to agree. */
QRectF FractalWidget::selection() const
{
    if (!m_selecting)
        return QRectF();
    QRectF r = QRectF(m_selectionStart, m_mousePosition).normalized();
    if (height() <= 0 || r.height() <= 0 || r.width() <= 0)
        return QRectF();
    const qreal want = (qreal)width() / (qreal)height();
    QPointF middle = r.center();
    qreal w = r.width(), h = r.height();
    if (w / h < want)
        w = h * want;
    else
        h = w / want;
    return QRectF(middle.x() - w / 2, middle.y() - h / 2, w, h);
}

void FractalWidget::mousePressEvent(QMouseEvent *event)
{
    m_mousePosition = event->pos();
    if (uih_selectionzoom_mode && event->button() == Qt::LeftButton) {
        m_selectionStart = event->pos();
        m_selecting = true;
        update();
        event->accept();
        return;
    }
    event->ignore();
}

void FractalWidget::mouseReleaseEvent(QMouseEvent *event)
{
    m_mousePosition = event->pos();
    if (m_selecting) {
        QRectF r = selection();
        m_selecting = false;
        update();
        if (m_uih != NULL && !r.isNull()) {
            /* Into device pixels, which is what the engine measures the image
             * in. */
            const qreal ratio = devicePixelRatioF();
            uih_zoomrectangle(m_uih, qRound(r.left() * ratio),
                              qRound(r.top() * ratio),
                              qRound(r.right() * ratio),
                              qRound(r.bottom() * ratio));
        }
        event->accept();
        return;
    }
    event->ignore();
}

void FractalWidget::mouseMoveEvent(QMouseEvent *event)
{
    m_mousePosition = event->pos();
    if (m_selecting) {
        update();
        event->accept();
        return;
    }
    event->ignore();
}

void FractalWidget::wheelEvent(QWheelEvent *event)
{
    m_mousePosition = event->position();
    event->ignore();
}
