#!/usr/bin/env roseus

(setq cnt 1000000)
(setq val 1.0)
(setq rate 5.3)

(format t "multiplication ~A times~%" cnt)
(setq s (unix::gettimeofday))
(dotimes (I cnt)
  (setq val (* val rate)))
(setq e (unix::gettimeofday))
(print val)
(format t ";; time ~A ms~%"
	(+ (* 1000.0 (- (car e) (car s)))
	   (/ (- (cadr e) (cadr s)) 1000.0)))

(format t "addition ~A times~%" cnt)
(setq val 1.0)
(setq s (unix::gettimeofday))
(dotimes (I cnt)
  (setq val (+ val rate)))
(setq e (unix::gettimeofday))
(print val)
(format t ";; time ~A ms~%"
	(+ (* 1000.0 (- (car e) (car s)))
	   (/ (- (cadr e) (cadr s)) 1000.0)))


(exit)
