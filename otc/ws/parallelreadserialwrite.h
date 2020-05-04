#ifndef PARALLEL_READ_SERIAL_WRITE_H
#define PARALLEL_READ_SERIAL_WRITE_H

#include <iostream>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <iomanip>
#include "assert.hh"

namespace otc {

class ParallelReadSerialWrite {
    public:
    // The following four variables are the core, they are protectedb by  mutex
    std::mutex mutex;
    std::size_t num_readers_working = 0;
    std::size_t num_writers_waiting = 0;
    void * the_working_writer = nullptr;
    // these conditional_variables keep us from having threads poll for their next available chance to run
    std::condition_variable no_readers_working_cond_var;
    std::condition_variable writer_released_cond_var;
    // a name for debugging, logging purposes
    const std::string state_name;

    explicit ParallelReadSerialWrite(const std::string & name)
        :state_name(name) {
        std::size_t num_readers_working = 0;
        std::size_t num_writers_waiting = 0;
        void * the_working_writer = nullptr;
    }
    bool write_possible() const {
        return num_readers_working == 0;
    }
    bool read_possible() const {
        return (num_writers_waiting == 0) && (the_working_writer == nullptr);
    }

    ParallelReadSerialWrite(const ParallelReadSerialWrite &) = delete;
    ParallelReadSerialWrite(const ParallelReadSerialWrite &&) = delete;
    ParallelReadSerialWrite & operator=(const ParallelReadSerialWrite &) = delete;
};

class ReadWriteMutex {
    protected:
    ParallelReadSerialWrite & shared;
    const char * const name; // just for logging and and debugging
    bool at_work = false;

    explicit ReadWriteMutex(ParallelReadSerialWrite & prsw,
                            const char * name_prefix)
        :shared(prsw),
        name(name_prefix) {
    }



    public:
    
    // Returns true if there are no writers waiting or working.
    // If update_bookkeeping is true, then it assumes that a response
    //    of true will cause the calling wrapper's thread to proceed.
    //    In this case the bookkeeping variables (num_readers_working)
    //    are updated to reflect the presence of another active reading thread.
    // If update_bookkeeping is false, it just returns the response to can_read
    //    as a const operation.
    bool can_read(bool update_bookkeeping) {
        bool writers_block;
        {
            std::unique_lock<std::mutex> book_lock(shared.mutex);
            //debug_book_lock_held("start can_read  " + std::to_string(update_bookkeeping), std::cerr);
            writers_block = (shared.num_writers_waiting > 0) || (shared.the_working_writer != nullptr) ;
            if (update_bookkeeping) {
                if (!writers_block) {
                    shared.num_readers_working += 1;
                }
            }
            //debug_book_lock_held("end   can_read  " + std::to_string(update_bookkeeping), std::cerr);
        }
        if (writers_block) {
            return false;
        }
        if (update_bookkeeping) {
            at_work = true;
        }
        return true;
    }
    
    // Must be called by every read-only wrapper that has obtained true from `can_read(true)`
    // This call decrements the num_readers_working count, and notifies 
    //    writers if that count has dropped to 0.
    void done_reading() {
        if (!at_work) { // no bookkeeping to be done if there was an excep
            return;
        }
        bool is_last_reader; //, writers_block;
        {
            std::unique_lock<std::mutex> book_lock(shared.mutex);
            //debug_book_lock_held("start done_read ", std::cerr);
            is_last_reader = (shared.num_readers_working == 1);
            //writers_block = (shared.num_writers_waiting > 0);
            assert(shared.the_working_writer == nullptr);
            shared.num_readers_working -= 1;
            //debug_book_lock_held("end   done_read ", std::cerr);
        }
        if (is_last_reader) {
            shared.no_readers_working_cond_var.notify_one();
        }
    }

    // for debugging, logging purposes
    std::string get_thread_wrapper_id() const {
        std::string x = name;
        return x;
    }
    ReadWriteMutex(const ReadWriteMutex &) = delete;
    ReadWriteMutex(const ReadWriteMutex &&) = delete;
    ReadWriteMutex operator=(const ReadWriteMutex &) = delete;
};

class  ReadMutexWrapper: public ReadWriteMutex {
    public:
    // the ctor blocks until a read-only operation can proceed
    explicit ReadMutexWrapper(ParallelReadSerialWrite & prsw)
        :ReadWriteMutex(prsw, "R") {
        while (!can_read(true)) {
            {
            }
            {
                std::unique_lock<std::mutex> ssul(shared.mutex);
                shared.writer_released_cond_var.wait(ssul, [&](){
                        return shared.read_possible(); //this_r_m_w->can_read(false);
                    });
            }
            {
            }
        }
    }
    ~ReadMutexWrapper() {
        done_reading();
    }
    ReadMutexWrapper(const ReadMutexWrapper &) = delete;
    ReadMutexWrapper(const ReadMutexWrapper &&) = delete;
    ReadMutexWrapper operator=(const ReadMutexWrapper &) = delete;
};

class  WriteMutexWrapper: public ReadWriteMutex {
    bool already_logged_as_waiting = false;
    public:
    bool can_write(bool update_bookkeeping) {
        bool readers_block, writers_block;
        {
            std::unique_lock<std::mutex> book_lock(shared.mutex);
            //debug_book_lock_held("start can_write " + std::to_string(update_bookkeeping) , std::cerr);
            readers_block = (shared.num_readers_working > 0) ;
            writers_block = false;
            if (readers_block) {
                if (update_bookkeeping && (!already_logged_as_waiting)) {
                    already_logged_as_waiting = true;
                    shared.num_writers_waiting += 1;
                }    
            } else {
                if (shared.the_working_writer != nullptr) {
                    if (update_bookkeeping && (!already_logged_as_waiting)) {
                        already_logged_as_waiting = true;
                        shared.num_writers_waiting += 1;
                    }
                    writers_block = true;
                } else {
                    if (update_bookkeeping) {
                        shared.the_working_writer = (void *)this;
                        if (already_logged_as_waiting) {
                            shared.num_writers_waiting -= 1;
                        }
                    }
                }
            }
            //debug_book_lock_held("end   can_write " + std::to_string(update_bookkeeping), std::cerr);
        }
        bool cw = !(readers_block || writers_block);
        if (cw && update_bookkeeping) {
            at_work = true;
        }
        return cw;
    }
    
    void done_writing() {
        if (!at_work) {
            return;
        }
        bool is_last_writer;
        {
            std::unique_lock<std::mutex> book_lock(shared.mutex);
            //debug_book_lock_held("start done_write", std::cerr);
            assert(shared.the_working_writer == (void *)this);
            shared.the_working_writer = nullptr;
            is_last_writer = (shared.num_writers_waiting == 0);
            //debug_book_lock_held("end   done_write", std::cerr);
        }
        if (is_last_writer) {
            {
            }
           shared.writer_released_cond_var.notify_all();
        } else {
            {
            }
            shared.no_readers_working_cond_var.notify_one();
        }
    }
    // the ctor blocks until a write operation can proceed
    explicit WriteMutexWrapper(ParallelReadSerialWrite & prsw)
        :ReadWriteMutex(prsw, "W"),
        already_logged_as_waiting(false) {
        while (!can_write(true)) { 
            {
            }
            {
                //WriteMutexWrapper * this_w_m_w = this;
                std::unique_lock<std::mutex> ssul(shared.mutex);
                shared.no_readers_working_cond_var.wait(ssul,
                                                        [&](){
                        return shared.write_possible();
                    });
            }
            {
            }
        }
    }
    ~WriteMutexWrapper() {
        done_writing();
    }
    
    WriteMutexWrapper(const WriteMutexWrapper &) = delete;
    WriteMutexWrapper(const WriteMutexWrapper &&) = delete;
    WriteMutexWrapper operator=(const WriteMutexWrapper &) = delete;

};


#if 0
void read_thread(ParallelReadSerialWrite * prsw, int tn) {
    try {
        int num_millisecond_sleep = 50;
        int n = 10;
        while (n >= 0) {
            std::this_thread::sleep_for(std::chrono::milliseconds(num_millisecond_sleep));
            {
                ReadMutexWrapper rmw(*prsw, tn);
                {
                    std::unique_lock<std::mutex> cout_lock(ParallelReadSerialWrite::cout_mutex);
                    std::cout << rmw.get_thread_wrapper_id() << " about to sleep n = " << n <<  std::endl;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(num_millisecond_sleep/10));
            }
            --n;
        }
    }  catch (...) {
        std::unique_lock<std::mutex> cout_lock(ParallelReadSerialWrite::cout_mutex);
        std::cout << "R" << tn << ": read_thread exiting with an exception " << std::endl;
        throw;
    }
    {
        std::unique_lock<std::mutex> cout_lock(ParallelReadSerialWrite::cout_mutex);
        std::cout << "R" << tn << ": read_thread exiting" << std::endl;
    }
}

void write_thread(ParallelReadSerialWrite * prsw, int tn) {
    try {
        int num_millisecond_sleep = 60;
        int n = 10;
        while (n >= 0) {
            std::this_thread::sleep_for(std::chrono::milliseconds(num_millisecond_sleep));
            {
                WriteMutexWrapper wmw(*prsw, tn);
                {
                    std::unique_lock<std::mutex> cout_lock(ParallelReadSerialWrite::cout_mutex);
                    std::cout << wmw.get_thread_wrapper_id() << " about to sleep n = " << n << std::endl;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(num_millisecond_sleep/10));
                {
                    std::unique_lock<std::mutex> cout_lock(ParallelReadSerialWrite::cout_mutex);
                    std::cout << wmw.get_thread_wrapper_id() << " slept about to release wrapper " << n << std::endl;
                }
                
            }
            --n;
        }
    } catch (...) {
        std::unique_lock<std::mutex> cout_lock(ParallelReadSerialWrite::cout_mutex);
        std::cout << "W" << tn << ": write_thread exiting with an exception " << std::endl;
        throw;
    }
    {
        std::unique_lock<std::mutex> cout_lock(ParallelReadSerialWrite::cout_mutex);
        std::cout << "W" << tn << ": write_thread exiting" << std::endl;
    }
}
#endif

} // namespace otc
 
#endif
