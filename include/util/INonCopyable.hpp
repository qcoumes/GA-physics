#ifndef GA_PHYSICS_INONCOPYABLE_HPP
#define GA_PHYSICS_INONCOPYABLE_HPP

namespace mastercraft::util {
    
    /**
     * Common interface for non copyable classes.
     */
    class INonCopyable {
        public:
            
            /**
             * Delete copy constructor.
             */
            INonCopyable(const INonCopyable &) = delete;
            
            /**
             * Delete assignment operator.
             */
            INonCopyable &operator=(const INonCopyable &) = delete;
        
        protected:
            INonCopyable() = default;
            
            ~INonCopyable() = default;
    };
}

#endif //GA_PHYSICS_INONCOPYABLE_HPP
