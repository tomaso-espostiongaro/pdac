      SUBROUTINE myflush( uni )
        INTEGER :: uni
#if defined __SP4
        CALL flush_( uni )
        ! flush( uni )
#endif
        RETURN
      END SUBROUTINE
