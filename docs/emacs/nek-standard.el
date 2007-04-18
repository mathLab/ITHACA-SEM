(provide 'nek-standard)

; c-no-initializer-indent
; 
(defun c-no-initializer-indent (langelem)
  ;; using in a statement block intro...
  (save-excursion
    (let* ((curpos (point))
   retval)
      (setq retval (if (search-forward "{" (c-point 'eol) t) 
       (* -1 c-basic-offset) c-basic-offset))
      (goto-char curpos)
      retval)))

; my-c-mode-common-hook
; 
(defun my-c-mode-common-hook ()
  ;; my customizations for all of c-mode and related modes
  (c-set-style "gnu")
  (set 'c-basic-offset 4)
  (set 'c-label-minimum-indentation 0)
  (c-set-offset 'substatement-open 0)
  (setq-default indent-tabs-mode nil);  
  (c-set-offset 'statement-cont 'c-no-initializer-indent)
  )

(add-hook 'c-mode-common-hook 'my-c-mode-common-hook)

