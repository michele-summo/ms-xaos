/* The formula bar: a QTextEdit whose document is laid out by NameWrapLayout
 * below rather than by Qt's own layout, so that a long formula goes on to the
 * next line between names. See formulaedit.h. */

#include <QtWidgets>
#include <cmath>

#include "formulaedit.h"
#include "formulahelp.h"

/* A document layout that ends each line where formula_wrap says.
 *
 * QAbstractTextDocumentLayout is Qt's way of laying a document out
 * differently, and QTextEdit works with any: it asks the layout for sizes and
 * hit tests, and finds the cursor through each block's own QTextLayout, which
 * is what this fills. Each line is given its exact number of characters with
 * QTextLine::setNumColumns, which is what lets the break fall between names
 * rather than where the Unicode line breaking rules would put it. */
class NameWrapLayout : public QAbstractTextDocumentLayout
{
  public:
    explicit NameWrapLayout(QTextDocument *document)
        : QAbstractTextDocumentLayout(document)
    {
    }

    void draw(QPainter *painter, const PaintContext &context) override;
    int hitTest(const QPointF &point,
                Qt::HitTestAccuracy accuracy) const override;
    int pageCount() const override { return 1; }
    QSizeF documentSize() const override { return m_size; }
    QRectF frameBoundingRect(QTextFrame *) const override
    {
        return QRectF(QPointF(0, 0), m_size);
    }
    QRectF blockBoundingRect(const QTextBlock &block) const override;

  protected:
    void documentChanged(int, int, int) override { relayout(); }

  private:
    void relayout();
    /* The height of a block's lines together. */
    static qreal linesHeight(const QTextLayout *layout);

    QSizeF m_size;
};

void NameWrapLayout::relayout()
{
    QTextDocument *doc = document();
    const qreal margin = doc->documentMargin();
    /* The width the bar gives, or none yet: then nothing is broken. */
    const qreal width = doc->textWidth();
    const qreal avail = width > 0 ? qMax<qreal>(1, width - 2 * margin) : 1e9;
    const QFontMetricsF metrics(doc->defaultFont());

    qreal y = margin, widest = 0;
    for (QTextBlock block = doc->begin(); block.isValid();
         block = block.next()) {
        const QString text = block.text();
        const std::vector<int> ends = formula_wrap(
            reinterpret_cast<const char16_t *>(text.utf16()), (int)text.size(),
            avail, [&text, &metrics](int start, int end) {
                return metrics.horizontalAdvance(text.mid(start, end - start));
            });

        /* Anywhere, which is what makes the count below exact: setNumColumns
         * ends a line only where the wrap mode allows one to end, and under
         * any other mode goes on to the next space -- which in a formula is
         * usually the end of it. */
        QTextLayout *layout = block.layout();
        QTextOption option = doc->defaultTextOption();
        option.setWrapMode(QTextOption::WrapAnywhere);
        layout->setTextOption(option);
        layout->beginLayout();
        qreal height = 0;
        int start = 0;
        for (int end : ends) {
            QTextLine line = layout->createLine();
            if (!line.isValid())
                break;
            /* setNumColumns counts glyphs, so a font that draws two letters
             * as one could hand a line a character more than it was given;
             * it is then given one fewer until it fits. */
            int columns = end - start;
            line.setNumColumns(columns);
            while (line.textLength() > end - start && columns > 1)
                line.setNumColumns(--columns);
            line.setPosition(QPointF(margin, height));
            height += line.height();
            widest = qMax(widest, line.naturalTextWidth());
            start = line.textStart() + line.textLength();
        }
        layout->endLayout();
        layout->setPosition(QPointF(0, y));
        y += height;
    }
    m_size = QSizeF(width > 0 ? width : widest + 2 * margin, y + margin);
    emit documentSizeChanged(m_size);
    emit update();
}

qreal NameWrapLayout::linesHeight(const QTextLayout *layout)
{
    qreal height = 0;
    for (int i = 0; i < layout->lineCount(); i++)
        height += layout->lineAt(i).height();
    return height;
}

QRectF NameWrapLayout::blockBoundingRect(const QTextBlock &block) const
{
    if (!block.isValid())
        return QRectF();
    const QTextLayout *layout = block.layout();
    return QRectF(layout->position(),
                  QSizeF(m_size.width(), linesHeight(layout)));
}

void NameWrapLayout::draw(QPainter *painter, const PaintContext &context)
{
    bool known = false;
    int cursorWidth = property("cursorWidth").toInt(&known);
    if (!known)
        cursorWidth = 1;
    painter->setPen(context.palette.color(QPalette::Text));
    for (QTextBlock block = document()->begin(); block.isValid();
         block = block.next()) {
        QTextLayout *layout = block.layout();
        const int from = block.position(), length = block.length();
        /* The selection, cut to this block: the bar draws it by the format
         * QTextEdit hands over, highlight colours and all. */
        QList<QTextLayout::FormatRange> selections;
        for (const Selection &s : context.selections) {
            const int a = qMax(s.cursor.selectionStart() - from, 0);
            const int b = qMin(s.cursor.selectionEnd() - from, length);
            if (b > a) {
                QTextLayout::FormatRange range;
                range.start = a;
                range.length = b - a;
                range.format = s.format;
                selections << range;
            }
        }
        layout->draw(painter, QPointF(0, 0), selections, context.clip);
        const int at = context.cursorPosition - from;
        if (context.cursorPosition >= 0 && at >= 0 && at < length)
            layout->drawCursor(painter, QPointF(0, 0), at, cursorWidth);
    }
}

int NameWrapLayout::hitTest(const QPointF &point,
                            Qt::HitTestAccuracy /*accuracy*/) const
{
    /* The line level with the point, or the nearest one above or below: a
     * click in the margin still lands in the text, as it does in any
     * editor. */
    QTextBlock block = document()->begin();
    for (; block.isValid(); block = block.next()) {
        const QTextLayout *layout = block.layout();
        const QPointF at = point - layout->position();
        const int lines = layout->lineCount();
        for (int i = 0; i < lines; i++) {
            const QTextLine line = layout->lineAt(i);
            const bool last = i == lines - 1 && !block.next().isValid();
            if (at.y() < line.y() + line.height() || last)
                return block.position() + line.xToCursor(at.x());
        }
    }
    return qMax(0, document()->characterCount() - 1);
}

FormulaEdit::FormulaEdit(QWidget *parent) : QTextEdit(parent)
{
    setAcceptRichText(false);
    setTabChangesFocus(true);
    setLineWrapMode(QTextEdit::WidgetWidth);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    /* The document takes the layout over and deletes the one it had, and
     * QTextEdit follows the change of its own accord. */
    document()->setDocumentLayout(new NameWrapLayout(document()));
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    /* A line more or less is a new height, which the layout holding the bar
     * has to be told of. */
    connect(document()->documentLayout(),
            &QAbstractTextDocumentLayout::documentSizeChanged, this,
            [this](const QSizeF &) { updateGeometry(); });
}

void FormulaEdit::setTallest(int pixels)
{
    m_tallest = pixels;
    updateGeometry();
}

void FormulaEdit::setPreferredWidth(int pixels)
{
    m_preferredWidth = pixels;
    updateGeometry();
}

int FormulaEdit::wantedHeight() const
{
    const QMargins margins = contentsMargins();
    const int frame = margins.top() + margins.bottom();
    const qreal margin = document()->documentMargin();
    const int line = (int)std::ceil(QFontMetricsF(document()->defaultFont())
                                        .lineSpacing() +
                                    2 * margin) +
                     frame;
    const int text =
        (int)std::ceil(document()->documentLayout()->documentSize().height()) +
        frame;
    return qBound(line, text, qMax(line, m_tallest));
}

QSize FormulaEdit::sizeHint() const
{
    return QSize(m_preferredWidth, wantedHeight());
}

QSize FormulaEdit::minimumSizeHint() const
{
    return QSize(100, wantedHeight());
}

void FormulaEdit::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_Return || event->key() == Qt::Key_Enter) {
        event->ignore();
        return;
    }
    QTextEdit::keyPressEvent(event);
}

void FormulaEdit::insertFromMimeData(const QMimeData *source)
{
    if (!source || !source->hasText())
        return;
    QString text = source->text();
    text.replace(QLatin1String("\r\n"), QLatin1String(" "));
    for (QChar &ch : text)
        if (ch == '\n' || ch == '\r' || ch == '\t' ||
            ch == QChar::ParagraphSeparator || ch == QChar::LineSeparator)
            ch = ' ';
    insertPlainText(text);
}

/* --- the formulas used before --------------------------------------------- */

static const char history_key[] = "Formulas/UserFormulas";
static const int history_most = 10;

QStringList formula_history()
{
    /* Once each, keeping the newest: a history written before the dialog took
     * this over can hold the same formula several times, one per position or
     * render that ran it. */
    QSettings settings;
    QStringList used = settings.value(history_key).toStringList();
    used.removeDuplicates();
    return used;
}

void formula_history_add(const QString &formula)
{
    if (formula.trimmed().isEmpty())
        return;
    /* Once each: a formula used again moves to the front rather than
     * appearing twice, which in a list of ten would crowd out another. */
    QStringList used = formula_history();
    used.removeAll(formula);
    used.prepend(formula);
    while (used.size() > history_most)
        used.removeLast();
    QSettings settings;
    settings.setValue(history_key, used);
}

/* The list in the history window, which words its entries again whenever its
 * width changes: how much of a long formula fits is a question of the width,
 * and Qt's own eliding says nothing about how long the formula was. */
class HistoryList : public QListWidget
{
  public:
    explicit HistoryList(QWidget *parent) : QListWidget(parent)
    {
        setTextElideMode(Qt::ElideNone);
        setWordWrap(false);
        setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        setUniformItemSizes(true);
    }

    void relabel()
    {
        const QFontMetricsF metrics(font());
        /* The item's own padding and the focus frame, on both sides. */
        const qreal avail =
            viewport()->width() - 4 * metrics.horizontalAdvance(QChar(' '));
        for (int i = 0; i < count(); i++) {
            QListWidgetItem *entry = item(i);
            const QString formula = entry->data(Qt::UserRole).toString();
            const QString suffix =
                QString("... ") +
                QObject::tr("(%1 characters)").arg(formula.size());
            const std::u16string shown = formula_history_label(
                formula.toStdU16String(), suffix.toStdU16String(), avail,
                [&metrics](const std::u16string &text) {
                    return metrics.horizontalAdvance(
                        QString::fromStdU16String(text));
                });
            entry->setText(QString::fromStdU16String(shown));
        }
    }

  protected:
    void resizeEvent(QResizeEvent *event) override
    {
        QListWidget::resizeEvent(event);
        relabel();
    }
};

void formula_history_open(QWidget *parent,
                          const std::function<void(const QString &)> &chosen)
{
    QDialog *window = new QDialog(parent);
    window->setAttribute(Qt::WA_DeleteOnClose);
    window->setWindowTitle(QObject::tr("Formulas used before"));

    const QStringList used = formula_history();
    QLabel *about = new QLabel(
        used.isEmpty()
            ? QObject::tr("No formula has been used yet.")
            : QObject::tr("Newest first. Double-click one, or pick it and "
                          "press Use, to put it in the formula bar."),
        window);
    about->setWordWrap(true);

    HistoryList *list = new HistoryList(window);
    for (const QString &formula : used) {
        QListWidgetItem *entry = new QListWidgetItem(formula, list);
        entry->setData(Qt::UserRole, formula);
    }

    QDialogButtonBox *buttons =
        new QDialogButtonBox(QDialogButtonBox::Cancel, window);
    QPushButton *use =
        buttons->addButton(QObject::tr("Use"), QDialogButtonBox::AcceptRole);
    /* Use takes the entry that is highlighted, which is the selected one: the
     * current one, with the focus frame, can be another -- after a Ctrl-click
     * that clears the selection, or a screen reader's select. */
    use->setEnabled(false);
    QObject::connect(list, &QListWidget::itemSelectionChanged, use,
                     [use, list]() {
                         use->setEnabled(!list->selectedItems().isEmpty());
                     });

    auto take = [window, chosen](QListWidgetItem *entry) {
        chosen(entry->data(Qt::UserRole).toString());
        window->accept();
    };
    QObject::connect(list, &QListWidget::itemDoubleClicked, window, take);
    QObject::connect(buttons, &QDialogButtonBox::accepted, window,
                     [take, list]() {
                         const QList<QListWidgetItem *> picked =
                             list->selectedItems();
                         if (!picked.isEmpty())
                             take(picked.first());
                     });
    QObject::connect(buttons, &QDialogButtonBox::rejected, window,
                     &QDialog::reject);

    QVBoxLayout *layout = new QVBoxLayout(window);
    layout->addWidget(about);
    layout->addWidget(list);
    layout->addWidget(buttons);

    /* Wide enough for a good part of a long formula, and within a 720p
     * screen like the dialog it opens from. */
    const QRect area = window->screen()->availableGeometry();
    window->resize(qMin(800, area.width() - 80), qMin(420, area.height() - 80));
    if (!used.isEmpty())
        list->setCurrentRow(0);
    window->open();
}
