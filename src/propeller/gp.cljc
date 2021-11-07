(ns propeller.gp
  (:require [clojure.string]
            [clojure.pprint]
            [propeller.genome :as genome]
            [propeller.variation :as variation]
            [propeller.push.instructions.bool]
            [propeller.push.instructions.character]
            [propeller.push.instructions.code]
            [propeller.push.instructions.input-output]
            [propeller.push.instructions.numeric]
            [propeller.push.instructions.polymorphic]
            [propeller.push.instructions.string]
            [propeller.push.instructions.vector]))

(defn report
  "Reports information each generation."
  [pop generation argmap]
  (let [best (first pop)]
    (clojure.pprint/pprint {:generation            generation
                            :best-plushy           (:plushy best)
                            :best-program          (genome/plushy->push (:plushy best) argmap)
                            :best-total-error      (:total-error best)
                            :best-errors           (:errors best)
                            :best-behaviors        (:behaviors best)
                            :genotypic-diversity   (float (/ (count (distinct (map :plushy pop))) (count pop)))
                            :behavioral-diversity  (float (/ (count (distinct (map :behaviors pop))) (count pop)))
                            :average-genome-length (float (/ (reduce + (map count (map :plushy pop))) (count pop)))
                            :average-total-error   (float (/ (reduce + (map :total-error pop)) (count pop)))})
    (println)))

(defn track-solutions
  "Track all solutions we've found so far.
   struct-solutions is a map where the keys are the hash values of solution
   programs, and the values are how many times we've seen them so far."
  [argmap error-function inputs population struct-solutions behavior-freq]
  (let [solutions (filter #(zero? (:total-error %)) population)]
    (loop [solutions solutions
           structs struct-solutions
           behavior-freq behavior-freq]
      (if (empty? solutions)
        [structs behavior-freq]
        (let [push-prog (genome/plushy->push (:plushy (first solutions)))
              prog-hash (hash push-prog)
              io-map {:inputs inputs
                      :outputs (repeat (count inputs) 0)}
              behavior (:behaviors (error-function argmap
                                                   io-map
                                                   (first solutions)))
              struct-map-for-behavior (get behavior-freq behavior {})]
          (recur (rest solutions)
                 (if (get structs prog-hash)
                   (update structs prog-hash inc)
                   (assoc structs prog-hash 1))
                 (assoc behavior-freq
                        behavior
                        (if (get struct-map-for-behavior prog-hash)
                          (update struct-map-for-behavior prog-hash inc)
                          (assoc struct-map-for-behavior prog-hash 1)))))))))

(defn end-report
  "Does calculations of L and b from paper"
  [struct-solutions behavior-freq]
  (println)
  (if (zero? (count behavior-freq))
    (do
      (println "No solutions")
      (println "L: 0")
      (println "n_struct: 0")
      (println "b: nil"))
    (let [modal-behavior (apply max-key
                                #(count (get behavior-freq %))
                                (keys behavior-freq))
          L (count (get behavior-freq modal-behavior))
          n_struct (count struct-solutions)
          b (- 2.0 (/ n_struct L))]
      (println "L:" L)
      (println "n_struct:" n_struct)
      (println "b:" b))))


(defn gp
  "Main GP loop."
  [{:keys [population-size max-generations error-function instructions
           max-initial-plushy-size solution-error-threshold random-inputs-for-generalizability mapper]
    :or   {solution-error-threshold 0.0
           ;; The `mapper` will perform a `map`-like operation to apply a function to every individual
           ;; in the population. The default is `map` but other options include `mapv`, or `pmap`.
           mapper #?(:clj pmap :cljs map)}
    :as   argmap}]
  ;;
  (prn {:starting-args (update (update argmap :error-function str) :instructions str)})
  (println)
  ;;
  (loop [generation 0
         population (repeatedly population-size
                                #(hash-map :plushy (genome/make-random-plushy
                                                    instructions
                                                    max-initial-plushy-size)))
         struct-solutions {}
         behavior-freq {}]
    (let [evaluated-pop (sort-by :total-error
                                 (mapper
                                   (partial error-function argmap (:training-data argmap))
                                   population))
          best-individual (first evaluated-pop)
          [struct-solutions behavior-freq] (track-solutions argmap
                                                            error-function
                                                            random-inputs-for-generalizability
                                                            evaluated-pop
                                                            struct-solutions
                                                            behavior-freq)]
      (if (:custom-report argmap)
        ((:custom-report argmap) evaluated-pop generation argmap)
        (report evaluated-pop generation argmap))
      (println "Struct Solutions Map:" (pr-str struct-solutions))
      (println "Behavioral Frequencies Map:" (pr-str behavior-freq))
      (println "Number of solutions:" (apply + (vals struct-solutions)))
      (println "Number of unique behaviors:" (count behavior-freq))
      (cond
        ;; Success on training cases is verified on testing cases
        ; (<= (:total-error best-individual) solution-error-threshold)
                false ;; This is to make it so we don't stop when finding a solution
        (do (prn {:success-generation generation})
            (prn {:total-test-error
                  (:total-error (error-function argmap (:testing-data argmap) best-individual))}))
        ;;
        (>= generation max-generations)
        (end-report struct-solutions behavior-freq)
        ;;
        :else (recur (inc generation)
                     (if (:elitism argmap)
                       (conj (repeatedly (dec population-size)
                                         #(variation/new-individual evaluated-pop argmap))
                             (first evaluated-pop))
                       (repeatedly population-size
                                   #(variation/new-individual evaluated-pop argmap)))
                     struct-solutions
                     behavior-freq)))))
