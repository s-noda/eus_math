#!/usr/bin/env roseus

(defvar *cv-plugin*
  (labels
   ((library_search
     (str &key
	  (depth 0)
	  colon-pos lib-path)
     (format t "  [~A] target=" depth)
     (cond
      ((eq (length str) 0)
       (format t "~% euscv_plugin_not_found in eus_cv.l~%")
       (exit -1))
      ((and (setq colon-pos (or (position #\: str) (length str)))
	    (setq lib-path (subseq str 0 colon-pos))
	    (setq lib-path
		  (if (eq (aref lib-path (- (length lib-path) 1)) #\/)
		      (subseq lib-path 0 (- (length lib-path) 1))
		    lib-path))
	    (probe-file (setq lib-path
			      (print (format nil "~A/libeus_cv.so" lib-path)))))
       (load-foreign lib-path))
      (t
       (library_search (subseq str (min (length str) (+ colon-pos 1)))
		       :depth (+ depth 1))))))
   (library_search (format nil "~A:~A/lib"
			   (unix:getenv "LD_LIBRARY_PATH")
			   (read-line (piped-fork "rospack find eus_math"))))))

(defforeign _cv-calc-centroid *cv-plugin* "cv_calc_centroid" (:integer :integer :string :string) :integer)
(defforeign _cv-calc-variance *cv-plugin* "cv_calc_variance" (:integer :integer :string :string :string) :integer)
(defforeign _cv-calc-gauss *cv-plugin* "cv_calc_gauss" (:integer :string :string :string) :float)
(defforeign _cv-calc-bic *cv-plugin* "cv_calc_bic" (:integer :integer :string :string :string) :float)
;;
(defforeign _cv-kmeans *cv-plugin* "cv_kmeans" (:integer :integer :integer :string :string :string :integer :double :integer) :integer)
(defforeign _cv-xmeans *cv-plugin* "cv_xmeans" (:integer :integer :integer :string :string :string :string :integer :double :integer :string :double :double :integer) :integer)
(defforeign _cv-hcluster *cv-plugin* "cv_hcluster" (:integer :integer :string :string :integer :integer :double :integer) :integer)

(defun cv-calc-centroid
  (p
   &key
   (d 1)
   (n (/ (length p) d))
   (points-str (dvector2float-bytestring p))
   (ret (instantiate string (* 4 d))))
  (_cv-calc-centroid n d points-str ret)
  (float-bytestring2dvector ret))

(defun cv-calc-variance
  (p
   &key
   (d 1)
   (n (/ (length p) d))
   (points-str (dvector2float-bytestring p))
   (c-str (let* ((buf (instantiate string (* 4 d))))
	    (_cv-calc-centroid n d points-str buf)
	    buf))
   (ret (instantiate string (* 4 d))))
  (_cv-calc-variance n d points-str c-str ret)
  (float-bytestring2dvector ret))

(defun cv-calc-gauss
  (p
   &key
   (points)
   (d (length p))
   (n (/ (length points) d))
   (points-str (dvector2float-bytestring points))
   (c-str (let* ((buf (instantiate string (* 4 d))))
	    (_cv-calc-centroid n d points-str buf)
	    buf))
   (var-str (let* ((buf (instantiate string (* 4 d))))
	      (_cv-calc-variance n d points-str c-str buf)
	      buf))
   )
  (_cv-calc-gauss d points-str c-str var-str))

(defun cv-calc-bic
  (p
   &key
   (d 1)
   (n (/ (length p) d))
   (points-str (dvector2float-bytestring p))
   (c-str (let* ((buf (instantiate string (* 4 d))))
	    (_cv-calc-centroid n d points-str buf)
	    buf))
   (var-str (let* ((buf (instantiate string (* 4 d))))
	      (_cv-calc-variance n d points-str c-str buf)
	      buf))
   )
  (_cv-calc-bic n d points-str c-str var-str))

(defun cv-xmeans
  (&key
   (sample-n 10)
   (length 3)
   (points (instantiate float-vector (* sample-n length)))
   (id (instantiate integer-vector sample-n))
   (centers (instantiate float-vector (* sample-n length)))
   (max-iter 10)
   (eps 0.01)
   ;;
   (points-str (dvector2float-bytestring points))
   (id-str (lvector2integer-bytestring id))
   (centers-str (dvector2float-bytestring centers))
   (var-str (copy-seq centers-str))
   (depth-str (instantiate string 4))
   ;;
   (max-depth 1000)
   (bic-stop-thre 0)
   (c-dist-thre 0)
   ;;
   (debug? t)
   )
  (_cv-xmeans 0 sample-n length points-str id-str centers-str var-str
	      max-iter eps max-depth depth-str bic-stop-thre c-dist-thre
	      (if debug? 1 0))
  ;;
  (if points (float-bytestring2dvector points-str points))
  (if centers (float-bytestring2dvector centers-str centers))
  (if id (integer-bytestring2lvector id-str id))
  ;;
  ;; (print points) (print centers) (print id)
  )

(defun cv-kmeans
  (&key
   (k 3)
   (sample-n 10)
   (length 3)
   (points (instantiate float-vector (* sample-n length)))
   (id (instantiate integer-vector sample-n))
   (centers (instantiate float-vector (* k length)))
   (max-iter 10)
   (eps 0.01)
   ;;
   (points-str (dvector2float-bytestring points))
   (id-str (lvector2integer-bytestring id))
   (centers-str (dvector2float-bytestring centers))
   ;;
   (debug? t)
   )
  (_cv-kmeans k sample-n length points-str id-str centers-str max-iter eps (if debug? 1 0))
  ;;
  (if points (float-bytestring2dvector points-str points))
  (if centers (float-bytestring2dvector centers-str centers))
  (if id (integer-bytestring2lvector id-str id))
  ;;
  ;; (print points) (print centers) (print id)
  )

(defun cv-hcluster
  (&key
   (sample-n 10)
   (length 3)
   (points (instantiate float-vector (* sample-n length)))
   (centers (instantiate float-vector (* sample-n length)))
   (branch 2)
   (max-iter 10)
   (eps 0.2)
   ;;
   (points-str (dvector2float-bytestring points))
   (centers-str (dvector2float-bytestring centers))
   ;;
   (debug? t)
   cnt
   )
  (setq
   cnt
   (_cv-hcluster sample-n length points-str centers-str
		 max-iter branch  eps (if debug? 1 0)))
  ;;
  (if points (float-bytestring2dvector points-str points))
  (if centers (float-bytestring2dvector centers-str centers))
  ;;
  ;; (print points) (print centers) (print id)
  cnt
  )

(defun gen-clusters
  (&key
   (sample-n 10)
   (length 3)
   (k 3)
   (ret (instantiate float-vector (* sample-n length)))
   p
   )
  (let* ((c 0))
    (dotimes (i k)
      (dotimes (j (ceiling (/ (* 1.0 sample-n) k)))
	(dotimes (h length)
	  (setq p (+ h (* j length)
		     (* i (ceiling (/ (* 1.0 sample-n) k)) length)))
	  (if (< p (length ret))
	      (setf (aref ret p) (+ c (random 1.0))))))
      (setq c (+ c 2))))
  ret)

(defun test-cv
  (&key
   (k 3)
   (sample-n 10)
   (length 3)
   (points)
   (id (instantiate integer-vector sample-n))
   (centers (instantiate float-vector (* k length)))
   (max-iter 10)
   (eps 0.01)
   ;;
   (test 'cv-kmeans)
   (rest-args nil)
   )
  (cond
   ((not points)
    (setq points (instantiate float-vector (* sample-n length)))
    (let* ((c 0))
      (dotimes (j length)
	(dotimes (i sample-n)
	  (setf (aref points (+ i (* j sample-n)))
		(+ c (random 1.0)))
	  )
	(setq c (+ c 100))))
    ))
  (cond
   ((eq test 'cv-kmeans)
    (cv-kmeans :k k :sample-n sample-n :length length :points points :id id :centers centers :max-iter max-iter :eps eps))
   ((eq test 'cv-hcluster)
    (cv-hcluster :sample-n sample-n :length length :points points :centers centers :max-iter max-iter :eps eps))
   ((eq test 'cv-xmeans)
    (setq centers (instantiate float-vector (* sample-n length)))
    (apply 'cv-xmeans :sample-n sample-n :length length :points points :id id :centers centers :max-iter max-iter :eps eps rest-args))
   )
  (print points)
  (print id)
  (print centers))

