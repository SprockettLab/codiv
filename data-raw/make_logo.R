suppressWarnings(suppressMessages({ library(ggplot2); library(hexSticker) }))

bg     <- "#f8fafc"; border <- "#2563eb"    # blue-600 (matches the morph theme)
tree   <- "#334155"; link   <- "#3b82f6"    # blue-500
wordcol<- "#0f172a"; boxfill <- "#3b82f6"

seg <- function(x1,y1,x2,y2) data.frame(x=x1,y=y1,xend=x2,yend=y2)

# RIGHT tree (symbiont): 3 tips, ((A,B),C), at y = 2,3,4 ; root stub reaches x=0.86
right <- do.call(rbind, list(
  seg(0.58,2,0.68,2), seg(0.58,3,0.68,3), seg(0.68,2,0.68,3),
  seg(0.68,2.5,0.76,2.5), seg(0.58,4,0.76,4), seg(0.76,2.5,0.76,4),
  seg(0.76,3.25,0.86,3.25)
))

# LEFT tree (host): 6 tips. Monophyletic matched clade at y 2,3,4 = ((2,3),4).
left <- do.call(rbind, list(
  seg(0.42,2,0.32,2), seg(0.42,3,0.32,3), seg(0.32,2,0.32,3),
  seg(0.32,2.5,0.24,2.5), seg(0.42,4,0.24,4), seg(0.24,2.5,0.24,4),
  seg(0.42,5,0.32,5), seg(0.42,6,0.32,6), seg(0.32,5,0.32,6),
  seg(0.24,3.25,0.16,3.25), seg(0.32,5.5,0.16,5.5), seg(0.16,3.25,0.16,5.5),
  seg(0.16,4.375,0.08,4.375), seg(0.42,1,0.08,1), seg(0.08,1,0.08,4.375),
  seg(0.02,2.6875,0.08,2.6875)
))
trees <- rbind(left, right)

links <- do.call(rbind, lapply(2:4, function(i) seg(0.42,i,0.58,i)))

# box covers the matched clade + right tree tips/cherry, but NOT the right root stub
boxes <- data.frame(xmin = 0.20, xmax = 0.80, ymin = 1.55, ymax = 4.45)

gg <- ggplot() +
  geom_rect(data = boxes, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill = boxfill, alpha = 0.22) +
  geom_segment(data = links, aes(x,y,xend=xend,yend=yend),
               color = link, linewidth = 0.8, lineend = "round") +
  geom_segment(data = trees, aes(x,y,xend=xend,yend=yend),
               color = tree, linewidth = 1.1, lineend = "round") +
  coord_cartesian(xlim = c(-0.05, 0.93), ylim = c(0.4, 6.6), expand = FALSE) +
  theme_void()

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
sticker(
  subplot = gg, s_x = 1, s_y = 1.17, s_width = 1.68, s_height = 1.16,
  package = "codiv", p_size = 26, p_x = 1, p_y = 0.50,
  p_color = wordcol, p_family = "sans", p_fontface = "bold",
  h_fill = bg, h_color = border, h_size = 1.5,
  dpi = 320, filename = "man/figures/logo.png"
)
cat("wrote man/figures/logo.png\n")
