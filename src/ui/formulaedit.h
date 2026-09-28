#ifndef FORMULAEDIT_H
#define FORMULAEDIT_H

#include <QStringList>
#include <QTextEdit>
#include <functional>

/* The formula bar of Fractal -> User formula: one formula, written on as many
 * lines as it takes.
 *
 * A formula is one line of text -- the parser reads no line breaks, and the
 * bar lets none in -- but one long enough is wider than any window, so it is
 * shown going on to the next line between names, never inside one (see
 * formula_wrap). Qt's own word wrap cannot do that: it breaks only where the
 * Unicode rules for prose allow, which in a formula is after some commas and
 * little else, and then cuts a name wherever it has to. So the bar lays its
 * text out itself.
 *
 * It is as tall as its lines, from one up to setTallest(), and scrolls past
 * that. Its size hint follows the text, so the dialog holding it grows and
 * shrinks with the formula.
 */
class FormulaEdit : public QTextEdit
{
  public:
    explicit FormulaEdit(QWidget *parent = nullptr);

    /* The tallest the bar may grow, in pixels, before its text scrolls. */
    void setTallest(int pixels);
    /* How wide the bar would like to be, which is what sizeHint gives. */
    void setPreferredWidth(int pixels);

    QSize sizeHint() const override;
    QSize minimumSizeHint() const override;

    /* The height the text needs, held between one line and the tallest. */
    int wantedHeight() const;

  protected:
    /* Return and Enter are the dialog's, which takes them as OK: a formula
     * has no second line to begin. */
    void keyPressEvent(QKeyEvent *event) override;
    /* Pasted or dropped text comes in as one line, its breaks and tabs made
     * spaces -- which the parser reads as nothing. */
    void insertFromMimeData(const QMimeData *source) override;

  private:
    int m_tallest = 400;
    int m_preferredWidth = 600;
};

/* The formulas used before: what Fractal -> User formula has been given, newest
 * first and each once, at most ten. Kept in the settings, and written only when
 * that dialog is accepted -- not when a position or a render runs usrform, which
 * is a formula being loaded rather than one being chosen. */
QStringList formula_history();
void formula_history_add(const QString &formula);

/* A window over parent listing them, each whole if it fits and otherwise
 * shortened with three dots and how many characters it has. The one picked
 * is handed to chosen. */
void formula_history_open(QWidget *parent,
                          const std::function<void(const QString &)> &chosen);

#endif // FORMULAEDIT_H
