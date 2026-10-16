#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "raylib.h"

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

#define WIDTH 800
#define HEIGHT 800
#define ANGLES 10
#define ENDTIME 1000

#define abs(a) ((a) > 0 ? (a) : -(a))
#define pixel(x, y, a) pixels[(y) * WIDTH * ANGLES + (x) * ANGLES + (a)]

double w, b, ener;
double t = 0;
double timestep = 0.1;
double max_r;
struct Pixel *pixels;

struct Pixel {
  double x, y, px, py;
  double endtime;
  int canal;
};


double ham(double x, double y, double px, double py)
{
  return (x*x + y*y + px*px + py*py)/2 - w*(x*py - y*px) + x*x*y - b*y*y*y;
}

void reset(float w_new, float b_new)
{
  w = w_new;
  b = b_new;
  max_r = (1 - w*w)*(1 - w*w)*(1 + b) + 0.1;

  ener = -(b+1)*(w*w - 1)*(w*w - 1)*(w*w - 1)/8;
  ener = ener * 1.13;

  double xmax = (1-w*w)*sqrt(2 + 3*b)/2 + 0.1;
  double xmin = -xmax;
  double ymax = (1-w*w)/(3*b) + 0.1;
  double ymin = (w*w-1)/2 - 0.1;

  for (int y = 0; y < HEIGHT; y++) {
    for (int x = 0; x < WIDTH; x++) {

      double xval = (x * (xmax - xmin)) / WIDTH + xmin;
      double yval = (y * (ymin - ymax)) / HEIGHT + ymax;

      double r = sqrt(2*ener + (w*w - 2*yval - 1)*xval*xval + (w*w + 2*b*yval - 1)*yval*yval);


      for (int i = 0; i < ANGLES; i++) {
        if (isnan(r)) {
          pixel(x, y, i).canal = -1;
          continue;
        } else {
          pixel(x, y, i).canal = 0;
        }
        double a = i*(2*M_PI / ANGLES);

        pixel(x, y, i).x = xval;
        pixel(x, y, i).y = yval;
        pixel(x, y, i).px = r*cos(a) - w*yval;
        pixel(x, y, i).py = r*sin(a) + w*xval;


        double H = ham(pixel(x, y, i).x, pixel(x, y, i).y, pixel(x, y, i).px, pixel(x, y, i).py);

        if (abs(H - ener) >= 1e-10) {
          printf("Fallo en %f, %f\n", xval, yval);
        }

        pixel(x, y, i).endtime = -1;
      }
    }
  }
}

static inline double
fx(double x, double y, double px, double py)
{
  return px + w*y;
}
static inline double
fy(double x, double y, double px, double py)
{
  return py - w*x;
}
static inline double
fpx(double x, double y, double px, double py)
{
  return py *w - x - 2*x*y;
}
static inline double
fpy(double x, double y, double px, double py)
{
  return -px*w - x*x - y + 3*b*y*y;
}

void step(struct Pixel *p)
{
  if (p->canal != 0) return;

  double x = p->x, y = p->y, px = p->px, py = p->py;

  double k1x = fx(x, y, px, py);
  double k1y = fy(x, y, px, py);
  double k1px = fpx(x, y, px, py);
  double k1py = fpy(x, y, px, py);

  double k2x = fx(x + k1x*timestep/2, y + k1y*timestep/2, px + k1px*timestep/2, py + k1py*timestep/2);
  double k2y = fy(x + k1x*timestep/2, y + k1y*timestep/2, px + k1px*timestep/2, py + k1py*timestep/2);
  double k2px = fpx(x + k1x*timestep/2, y + k1y*timestep/2, px + k1px*timestep/2, py + k1py*timestep/2);
  double k2py = fpy(x + k1x*timestep/2, y + k1y*timestep/2, px + k1px*timestep/2, py + k1py*timestep/2);

  double k3x = fx(x + k2x*timestep/2, y + k2y*timestep/2, px + k2px*timestep/2, py + k2py*timestep/2);
  double k3y = fy(x + k2x*timestep/2, y + k2y*timestep/2, px + k2px*timestep/2, py + k2py*timestep/2);
  double k3px = fpx(x + k2x*timestep/2, y + k2y*timestep/2, px + k2px*timestep/2, py + k2py*timestep/2);
  double k3py = fpy(x + k2x*timestep/2, y + k2y*timestep/2, px + k2px*timestep/2, py + k2py*timestep/2);

  double k4x = fx(x + k3x*timestep, y + k3y*timestep, px + k3px*timestep, py + k3py*timestep);
  double k4y = fy(x + k3x*timestep, y + k3y*timestep, px + k3px*timestep, py + k3py*timestep);
  double k4px = fpx(x + k3x*timestep, y + k3y*timestep, px + k3px*timestep, py + k3py*timestep);
  double k4py = fpy(x + k3x*timestep, y + k3y*timestep, px + k3px*timestep, py + k3py*timestep);

  p->x = x + (k1x + 2*k2x + 2*k3x + k4x) * timestep/ 6;
  p->y = y + (k1y + 2*k2y + 2*k3y + k4y) * timestep/ 6;
  p->px = px + (k1px + 2*k2px + 2*k3px + k4px) * timestep/ 6;
  p->py = py + (k1py + 2*k2py + 2*k3py + k4py) * timestep/ 6;

  double H = ham(p->x, p->y, p->px, p->py);

  if (abs(H - ener) >= 1e-4) {
    printf("Fallo: %f\n", H-ener);
    printf("Expected: %f\n", ener);
    printf("Got:      %f\n", H);
    printf("Position: (%f, %f, %f, %f)\n", p->x, p->y, p->px, p->py);
    printf("R:        %f\n", sqrt(p->x*p->x + p->y*p->y));
  }
  if (p->x*p->x + p->y*p->y >= max_r) {
    if (p->y > 0) p->canal = 3;
    else if (p->x > 0) p->canal = 2;
    else p->canal = 1;
  }
}

static inline
Color compute_color(int x, int y)
{
  if (pixel(x, y, 0).canal < 0) return (Color){ .r = 255, .g = 255, .b = 255, .a = 255 };

  int canal1 = 0;
  int canal2 = 0;
  int canal3 = 0;

  for (int i = 0; i < ANGLES; i++) {
    if (pixel(x, y, i).canal == 1) canal1++;
    if (pixel(x, y, i).canal == 2) canal2++;
    if (pixel(x, y, i).canal == 3) canal3++;
  }

  return (Color){
    .r = 255*canal1/ANGLES,
    .g = 255*canal2/ANGLES,
    .b = 255*canal3/ANGLES,
    .a = 255,
  };
}

void output(FILE *file)
{
  fprintf(file, "P3\n");
  fprintf(file, "\n");
  fprintf(file, "%d %d\n", WIDTH, HEIGHT);
  fprintf(file, "255\n");

  for (int y = 0; y < HEIGHT; y++) {
    for (int x = 0; x < WIDTH; x++) {
      Color c = compute_color(x, y);
      fprintf(file, "%d %d %d\n", c.r, c.g, c.b);
    }
  }
}

int main() {
  pixels = malloc(ANGLES*WIDTH*HEIGHT*sizeof(struct Pixel));

  char text[1024];

  InitWindow(WIDTH, HEIGHT, "Integration");

  Image img = {
    .data = malloc(WIDTH * HEIGHT * 4),
    .width = WIDTH,
    .height = HEIGHT,
    .format = PIXELFORMAT_UNCOMPRESSED_R8G8B8A8,
    .mipmaps = 1
  };

  if (!IsImageValid(img)) {
    fprintf(stderr, "Failed to create image\n");
    exit(1);
  }
  Texture2D texture = LoadTextureFromImage(img);

  SetTargetFPS(60);

  reset(0.1, 1./3);

  while (!WindowShouldClose()) {
    // STEP
    #pragma omp parallel for collapse(3) schedule(dynamic, 1)
    for (int y = 0; y < HEIGHT; y++) {
      for (int x = 0; x < WIDTH; x++) {
        for (int i = 0; i < ANGLES; i++) {
          step(&pixel(x, y, i));
        }
      }
    }
    t += timestep;

    // UPDATE IMAGE
    for (int y = 0; y < HEIGHT; y++) {
      for (int x = 0; x < WIDTH; x++) {
        Color *c = ((Color *)img.data) + y*WIDTH + x;
        *c = compute_color(x, y);
      }
    }
    UpdateTexture(texture, img.data);

    // DRAWING
    BeginDrawing();
      ClearBackground(RAYWHITE);
      DrawTexture(texture, 0, 0, WHITE);

      sprintf(text, "t: %f", t);
      DrawText(text, 0, 0, 20, BLACK);
      DrawFPS(0, 20);

      float w_f = w, b_f = b;
      GuiSliderBar((Rectangle){ 650, 10, 120, 20 }, "w", TextFormat("%.2f", w), &w_f, 0, 1);
      GuiSliderBar((Rectangle){ 650, 40, 120, 20 }, "b", TextFormat("%.2f", b), &b_f, 0, 1);
      if (w_f != w || b_f != b) {
        reset(w_f, b_f);
        t = 0;
      }

      if (GuiButton((Rectangle){650, 70, 120, 20 }, "Output file")) {
        char filename[1024];
        sprintf(filename, "generated/w%fbeta%f.ppm", w, 3*b);
        FILE *file = fopen(filename, "wb");
        output(file);
        fclose(file);
      }
    EndDrawing();
  }

  CloseWindow();

  // TODO(ernesto): Is it bad if we leak?
}
