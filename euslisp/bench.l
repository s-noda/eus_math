#!/usr/bin/env roseus

(setq cnt 100000)
(setq val 1.0)
(setq rate 5.3)
;; (setq ll (make-list cnt))

(defun dotimes2 (cnt val rate)
  (if (<= cnt 0) val (dotimes2 (- cnt 1) (* val rate) rate)))

(format t "multiplication ~A times~%" cnt)
(setq s (unix::gettimeofday))
(dotimes (i cnt)
;; (setq i 0) (while (< (setq i (+ i 1)) cnt)
;; (dolist (l ll)
  (setq val (* val rate)))
;; (dotimes2 cnt val rate)
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
