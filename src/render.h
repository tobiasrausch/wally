#ifndef RENDER_H
#define RENDER_H

#include <string>
#include <cstdint>
#include <map>

#include <blend2d/blend2d.h>

#include "font.h"

namespace wallysworld
{

  #ifndef WALLY_FONT_PX
  #define WALLY_FONT_PX 12.5f
  #endif

  #ifndef WALLY_FONT_REFSCALE
  #define WALLY_FONT_REFSCALE 0.4
  #endif

  inline BLFontFace&
  wallyFontFace() {
    static BLFontFace face = []() {
      BLFontData fd;
      fd.create_from_data(dejavuSansTTF, dejavuSansTTFLen);
      BLFontFace f;
      f.create_from_data(fd, 0);
      return f;
    }();
    return face;
  }

  inline const BLFont&
  wallyFont(float const px = WALLY_FONT_PX) {
    static std::map<int32_t, BLFont> cache;
    int32_t key = (int32_t) (px * 10.0f + 0.5f);
    std::map<int32_t, BLFont>::iterator it = cache.find(key);
    if (it != cache.end()) return it->second;
    BLFont f;
    f.create_from_face(wallyFontFace(), px);
    return cache.emplace(key, std::move(f)).first->second;
  }

  inline float
  scaleToPx(double const scale) {
    return (float) (scale * (WALLY_FONT_PX / WALLY_FONT_REFSCALE));
  }

  struct TextSize {
    int32_t width;
    int32_t height;
  };

  inline TextSize
  getTextSize(std::string const& text, float const px = WALLY_FONT_PX) {
    const BLFont& font = wallyFont(px);
    BLGlyphBuffer gb;
    gb.set_utf8_text(text.c_str());
    font.shape(gb);
    BLTextMetrics tm;
    font.get_text_metrics(gb, tm);
    TextSize ts;
    ts.width = (int32_t) (tm.advance.x + 0.5);
    ts.height = (int32_t) (px * 0.72f + 0.5f);
    return ts;
  }

  inline void
  drawText(BLContext& img, int32_t const x, int32_t const y, std::string const& text, BLRgba32 const& clr, float const px = WALLY_FONT_PX) {
    img.fill_utf8_text(BLPoint((double) x, (double) y), wallyFont(px), text.c_str(), SIZE_MAX, clr);
  }

  inline void
  drawTextBold(BLContext& img, int32_t const x, int32_t const y, std::string const& text, BLRgba32 const& clr, float const px = WALLY_FONT_PX) {
    const BLFont& font = wallyFont(px);
    img.fill_utf8_text(BLPoint((double) x, (double) y), font, text.c_str(), SIZE_MAX, clr);
    img.fill_utf8_text(BLPoint((double) x + 0.55, (double) y), font, text.c_str(), SIZE_MAX, clr);
  }

  inline void
  drawChip(BLContext& img, int32_t const x, int32_t const baseline, std::string const& text, BLRgba32 const& textClr, BLRgba32 const& bgClr, bool const bold = false) {
    TextSize ts = getTextSize(text);
    int32_t padX = 3;
    int32_t padY = 1;
    int32_t topY = baseline - ts.height - padY;
    int32_t w = ts.width + 2 * padX;
    int32_t h = ts.height + 2 * padY;
    img.fill_round_rect(BLRoundRect((double) x, (double) topY, (double) w, (double) h, 2.5, 2.5), bgClr);
    if (bold) drawTextBold(img, x + padX, baseline, text, textClr);
    else drawText(img, x + padX, baseline, text, textClr);
  }

  inline void
  drawLine(BLContext& img, double const x0, double const y0, double const x1, double const y1, BLRgba32 const& clr, double const w = 1.0) {
    img.set_stroke_width(w);
    img.stroke_line(BLPoint(x0 + 0.5, y0 + 0.5), BLPoint(x1 + 0.5, y1 + 0.5), clr);
  }

  inline void
  drawTextRotated(BLContext& img, double const ax, double const ay, std::string const& text, BLRgba32 const& clr, float const px = WALLY_FONT_PX) {
    img.save();
    img.translate(ax, ay);
    img.rotate(1.5707963267948966);  // +90 degrees
    img.fill_utf8_text(BLPoint(0.0, 0.0), wallyFont(px), text.c_str(), SIZE_MAX, clr);
    img.restore();
  }

}

#endif
