#lang racket

(require racket/draw
         plot/utils ; vector fns
         racket/generic)

(define WIDTH 320)
(define HEIGHT 240)


(define (mk-cam x y z) (vnormalize (vector x y z)))
(define cam (mk-cam 1 1 1))
(define right (vnormalize (vcross #(0 -1 0) cam)))
(define up (vnormalize (vcross right cam)))

;; these constants are used when calculating rays
;; -- needed to maintain proper aspect ratio
;;    bc both x and yoffset are multiplied by unit vector
;; ie if WIDTH > HEIGHT, then 0 < xoffset < 1
;;                        and 0 < yoffset < HEIGHT/WIDTH
(define DIV (if (> WIDTH HEIGHT) WIDTH HEIGHT))
(define (mk-ray x y)
  (define xoffset (/ x DIV))
  (define yoffset (/ y DIV))
  (vnormalize (v+ cam (v+ (v* right xoffset)
                          (v* up yoffset)))))

(define rays
  (for*/list ([x (in-range WIDTH)] [y (in-range HEIGHT)])
    (mk-ray x y)))
;    (define xoffset (/ x DIV))
;    (define yoffset (/ y DIV))
;    (vnormalize (v+ cam (v+ (v* right xoffset)
;                            (v* up yoffset))))))

(define-generics shape
  (intersect shape ray)
  (get-norm shape ray))

(struct Sphere (center radius)
  #:methods gen:shape
  [(define (intersect s ray)
     (define center (Sphere-center s))
     (define radius (Sphere-radius s))
     ;; proj is the projection of center onto ray (forming rt triangle)
     (define proj (vdot center ray))
     ;; neg proj means angle between vecs > 90, ie going away from each other
     (if (< proj 0) #f
                ;; r^2 + proj^2 = hypotenuse^2 of a rt triangle
         (let* ([hyp^2 (+ (* radius radius)
                          (* proj proj))]
                ;; (vdot cent cent) = (dist from 0 to center)^2
                [tmp (- hyp^2 (vdot center center))])
           ;; if this hyp^2 > (dist from 0 to center)^2, then the ray intersects
           ;; ow ray is outside the circle
           (if (< tmp 0) #f
               ;; writing out both sets of pythag eq and subtracting:
               ;; hyp^2 - (vdot cen cen) = tmp = r^2 - opp^2 = adj^2
               ;; where adj and opp are sides of a rt triangle where the hyp is
               ;; the center and the ray intersect point
               ;; so (sqrt tmp) gives adj, which when subtracted from proj gives
               ;; the dist to the intersection point
               ;; (it's easier to see if you draw it out on paper)
               (- proj (sqrt tmp))))))
   (define (get-norm s pt)
     (vnormalize (v- pt (Sphere-center s))))])
               
(define sph (Sphere #(4 12 10) 3))

(define lght #(10 0 0))

(define scene (make-bitmap WIDTH HEIGHT))

(define pix (make-bytes 32))

(send scene get-argb-pixels 0 0 2 2 pix)

pix 

(define red-pix (bytes 255 255 0 0))
(define green-pix (bytes 255 0 255 0))
(define blue-pix (bytes 128 0 0 255))
(define black-pix (bytes 255 0 0 0))

(for* ([x (in-range WIDTH)]
       [y (in-range HEIGHT)])
  (define ray (mk-ray x y))
  (define isect (intersect sph ray))
  (if isect
      (let ([isectpt (v* ray isect)])
        (define normray (get-norm sph isectpt))
        (define lightreflectray (vnormalize (v- lght isectpt)))
        (define intensity (abs (vdot lightreflectray normray)))
        (if (> intensity 1)
            (send scene set-argb-pixels x y 1 1 black-pix)
            (send scene set-argb-pixels x y 1 1 
                  (bytes (inexact->exact (round (* 255 intensity)))
                         (inexact->exact (round (* 255 intensity)))
                         0 0))))
      (send scene set-argb-pixels x y 1 1 black-pix)))
  
(send scene save-file "a.png" 'png)