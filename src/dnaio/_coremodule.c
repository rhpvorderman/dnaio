#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"         // PyMemberDef
#if defined(USE_SSE2)
    #include "ascii_check_sse2.h"
#else
    #include "ascii_check.h"
#endif
#include "_conversions.h"

static inline int 
record_ids_match_partial(
    char * header1, char * header2,
    size_t id1_length, size_t header2_length,
    int id1_ends_with_number) 
{
    if (header2_length < id1_length) {
        return 0;
    }
    char end = header2[id1_length];
    if ((end != 0) && (end != ' ') && (end != '\t')) {
        return 0;
    }
    char header2lastchar = header2[id1_length - 1];
    int id2_ends_with_number = (header2lastchar >= '1') && (header2lastchar <= '3');
    if (id1_ends_with_number && id2_ends_with_number) {
        id1_length -= 1;
    }
    return (memcmp(header1, header2, id1_length) == 0);
}

static inline int 
record_ids_match(char *header1, char *header2,
                 size_t header1_length) 
{
    size_t id2_length = strcspn(header2, " \t");
    char id2_end = header2[id2_length - 1];
    int id2_ends_with_number = (id2_end >= '1') && (id2_end <= '3');
    return record_ids_match_partial(header2, header1, id2_length, 
                                    header1_length, id2_ends_with_number);
}

typedef struct {
    PyObject_HEAD
    PyObject * name;
    PyObject * sequence;
    PyObject * qualities;
} SequenceRecord;

static void 
SequenceRecord_dealloc(SequenceRecord *self) {
    Py_CLEAR(self->name);
    Py_CLEAR(self->sequence);
    Py_CLEAR(self->qualities);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyTypeObject SequenceRecord_Type;

// NOTE: This function steals references. This is also how CPython tuples work.
// This means that when we create new name, sequence and qualities objects which
// will have their refcount automatically set to 1, we don't need to worry
// about incref or decref. There will be only one reference, and that will be
// owned by the sequence record object.
static inline PyObject *
new_sequence_record(PyObject *name, PyObject *sequence, PyObject *qualities){
    SequenceRecord *new_obj = PyObject_New(SequenceRecord, &SequenceRecord_Type);
    if (new_obj == NULL)
        return PyErr_NoMemory();
    new_obj->name = name;
    new_obj->sequence = sequence;
    new_obj->qualities = qualities;
    return (PyObject *)new_obj;
}

PyDoc_STRVAR(SequenceRecord__init____doc__,
"SequenceRecord(name, sequence, qualities = None)\n"
"--\n"
"\n"
"A named sequence with optional quality values.\n"
"This typically represents a record from a FASTA or FASTQ file.\n"
"The readers returned by `dnaio.open` yield objects of this type\n"
"when mode is set to ``\"r\"``\n"
"\n"
"Attributes:"
"    name (str): The read header\n"
"    sequence (str): The nucleotide (or amino acid) sequence\n"
"    qualities (str): None if no quality values are available\n"
"        (such as when the record comes from a FASTA file).\n"
"        If quality values are available, this is a string\n"
"        that contains the Phred-scaled qualities encoded as\n"
"        ASCII(qual+33) (as in FASTQ).\n"
"\n"
);

static PyObject *
SequenceRecord__new__(PyTypeObject *tp, PyObject *args, PyObject *kwargs) 
{
    PyObject *name = NULL;
    PyObject *sequence = NULL;
    PyObject *qualities = NULL;
    static char * _keywords[] = {"name", "sequence", "qualities", NULL};
    static char * _format = "O!O!|O:SequenceRecord";
    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, _format, _keywords,
        (PyObject *)&PyUnicode_Type, &name,
        (PyObject *)&PyUnicode_Type, &sequence,
        &qualities)) {
        return NULL;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(name)) {
        PyErr_SetString(PyExc_ValueError, "name must be a valid ASCII-string.");
        return NULL;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(sequence)) {
        PyErr_SetString(PyExc_ValueError, "sequence must be a valid ASCII-string.");
        return NULL;
    }
    if (qualities == Py_None) {
        qualities = NULL;
    }
    if (qualities != NULL) {
        if (!PyUnicode_CheckExact(qualities)) {
            PyErr_Format(
                PyExc_TypeError, 
                "qualities must be of type str, got: %s", 
                Py_TYPE(qualities)->tp_name);
            return NULL;
        }
        if (!PyUnicode_IS_COMPACT_ASCII(qualities)) {
            PyErr_SetString(PyExc_ValueError, 
                            "qualities must be a valid ASCII-string.");
            return NULL;
        }
        // Type already checked, can use unsafe macros here.
        if(PyUnicode_GET_LENGTH(sequence) != PyUnicode_GET_LENGTH(qualities)) {
            PyErr_Format(PyExc_ValueError,
                "Size of sequence and qualities do not match: %ld != %ld",
                PyUnicode_GET_LENGTH(sequence), PyUnicode_GET_LENGTH(qualities));
            return NULL;
        }
        Py_INCREF(qualities);
    }
    Py_INCREF(name);
    Py_INCREF(sequence);
    SequenceRecord * self = PyObject_New(SequenceRecord, tp);
    self->name = name;
    self->sequence=sequence;
    self->qualities=qualities;
    return (PyObject *)self;
};

// GETTERS AND SETTERS
static PyObject *
SequenceRecord_get_name(SequenceRecord *self, void *closure) 
{
    Py_INCREF(self->name);
    return self->name;
}

static PyObject *
SequenceRecord_get_sequence(SequenceRecord *self, void *closure) 
{
    Py_INCREF(self->sequence);
    return self->sequence;
}

static PyObject *
SequenceRecord_get_qualities(SequenceRecord *self, void *closure)
{
    if (self->qualities == NULL) {
        Py_RETURN_NONE;
    }
    Py_INCREF(self->qualities);
    return self->qualities;
}

static int
SequenceRecord_set_name(SequenceRecord *self, PyObject *value, void *closure)
{
    if (value == NULL){
        PyErr_SetString(PyExc_AttributeError, "name attribute cannot be deleted.");
        return -1;
    }
    if (!PyUnicode_CheckExact(value)) {
        PyErr_Format(PyExc_TypeError, 
                     "name must be of type str. Got %s", 
                     Py_TYPE(value)->tp_name);
        return -1;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(value)) {
        PyErr_SetString(PyExc_ValueError, "name must be a valid ASCII-string.");
        return -1;
    }
    PyObject * tmp = self->name;
    Py_INCREF(value);
    self->name = value;
    Py_DECREF(tmp);
    return 0;
}

static int
SequenceRecord_set_sequence(SequenceRecord *self, PyObject *value, void *closure)
{
    if (value == NULL){
        PyErr_SetString(PyExc_AttributeError, "sequence attribute cannot be deleted.");
        return -1;
    }
    if (!PyUnicode_CheckExact(value)) {
        PyErr_Format(PyExc_TypeError, 
                     "sequence must be of type str. Got %s", 
                     Py_TYPE(value)->tp_name);
        return -1;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(value)) {
        PyErr_SetString(PyExc_ValueError, "sequence must be a valid ASCII-string.");
        return -1;
    }
    PyObject * tmp = self->sequence;
    Py_INCREF(value);
    self->sequence = value;
    Py_DECREF(tmp);
    return 0;
}

static int
SequenceRecord_set_qualities(SequenceRecord *self, PyObject *value, void *closure)
{
    PyObject * tmp;
    if (value == NULL){
        PyErr_SetString(PyExc_AttributeError, "qualities attribute cannot be deleted.");
        return -1;
    }
    if (value == Py_None) {
        tmp = self->qualities;
        self->qualities = NULL;
        Py_DECREF(tmp);
        return 0;
    }
    if (!PyUnicode_CheckExact(value)) {
        PyErr_Format(PyExc_TypeError, 
                     "sequence must be of type str. Got %s", 
                     Py_TYPE(value)->tp_name);
        return -1;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(value)) {
        PyErr_SetString(PyExc_ValueError, "qualities must be a valid ASCII-string.");
        return -1;
    }
    tmp = self->qualities;
    Py_INCREF(value);
    self->qualities = value;
    Py_DECREF(tmp);
    return 0;
}

static PyGetSetDef SequenceRecord_properties[] = {
    {"name", (getter)SequenceRecord_get_name, (setter)SequenceRecord_set_name},
    {"sequence", (getter)SequenceRecord_get_sequence, (setter)SequenceRecord_set_sequence},
    {"qualities", (getter)SequenceRecord_get_qualities, (setter)SequenceRecord_set_qualities},
    {NULL}
};

// METHODS

static PyObject * 
SequenceRecord__repr__(SequenceRecord * self){
    const char * type_name = Py_TYPE(self)->tp_name;
    // Strip off module name from type name.
    char * type_name_after_dot = strchr(type_name, '.') + 1;
    if (self->qualities == NULL) {
        return PyUnicode_FromFormat("%s(%R, %R)", 
            type_name_after_dot, self->name, self->sequence);
    }
    return PyUnicode_FromFormat("%s(%R, %R, %R)", 
        type_name_after_dot, self->name, self->sequence, self->qualities);
}

static int 
SequenceRecord_equals(SequenceRecord * self, SequenceRecord * other)
{
    if (self->qualities == NULL) {
        if (other->qualities != NULL) {
            return 0;
        }
    return (PyObject_RichCompareBool(self->name, other->name, Py_EQ) && 
            PyObject_RichCompareBool(self->sequence, other->sequence, Py_EQ));
    }
    return (PyObject_RichCompareBool(self->name, other->name, Py_EQ) && 
            PyObject_RichCompareBool(self->sequence, other->sequence, Py_EQ) &&
            PyObject_RichCompareBool(self->qualities, other->qualities, Py_EQ)
            );
}
static PyObject *
SequenceRecord__richcompare__(SequenceRecord *self, SequenceRecord *other, int op)
{
    // This function is extremely generic to allow subtyping, reuse etc.
    if(Py_TYPE(self) != Py_TYPE(other)) {
        PyErr_Format(PyExc_TypeError, 
            "Can only compare objects of %R to objects of the same type. Got: %R.",
            Py_TYPE(self), Py_TYPE(other));
        return NULL;
    }
    if (op == Py_EQ){
        return PyBool_FromLong(SequenceRecord_equals(self, other));
    } else if (op == Py_NE) {
        return PyBool_FromLong(!SequenceRecord_equals(self, other));
    }
    else {
        return Py_NotImplemented;
    }
}

PyDoc_STRVAR(SequenceRecord_fastq_bytes__doc__,
"Return the entire FASTQ record as bytes which can be written\n"
"into a file.");

#define SEQUENCE_FASTQ_BYTES_METHODDEF    \
    {"fastq_bytes", (PyCFunction)(void(*)(void))SequenceRecord_fastq_bytes, \
     METH_FASTCALL | METH_KEYWORDS, SequenceRecord_fastq_bytes__doc__}

static PyObject *
SequenceRecord_fastq_bytes(SequenceRecord *self, 
                           PyObject *const *args, 
                           Py_ssize_t nargs, 
                           PyObject *kwnames){
    int two_headers = 0;
    if (nargs || kwnames) {
        if (nargs > 1) {
            PyErr_Format(PyExc_TypeError, 
                         "fastq_bytes() takes at most 1 argument (%ld given)", 
                          nargs);
            return NULL;
        }
        two_headers = PyObject_IsTrue(args[0]);
        if (kwnames) {
            Py_ssize_t nkwargs = PyTuple_GET_SIZE(kwnames);
            if ((nargs + nkwargs) > 1) {
                PyErr_Format(PyExc_TypeError, 
                         "fastq_bytes() takes at most 1 argument (%ld given)", 
                          nargs + nkwargs);
                return NULL;
            }
            PyObject * argname = PyTuple_GET_ITEM(kwnames, 0);
            if (strcmp(PyUnicode_DATA(argname), "two_headers") != 0) {
                PyErr_Format(PyExc_TypeError, "fastq_bytes() got an unexpected keyword argument %R", argname);
                return NULL;
            }
        }
    }
    if (self->qualities == NULL) {
        PyErr_SetString(PyExc_ValueError, 
        "Cannot create FASTQ bytes from a sequence without qualities.");
        return NULL;
    }
    Py_ssize_t name_length = PyUnicode_GET_LENGTH(self->name);
    Py_ssize_t sequence_length = PyUnicode_GET_LENGTH(self->sequence);
    Py_ssize_t qualities_length = PyUnicode_GET_LENGTH(self->qualities);
   
    char * name = PyUnicode_DATA(self->name);
    char * sequence = PyUnicode_DATA(self->sequence);
    char * qualities = PyUnicode_DATA(self->qualities);

    // Total size is name + sequence + qualities + 4 newlines + '+' and an
    // '@' to be put in front of the name.
    Py_ssize_t total_size = name_length + sequence_length + qualities_length + 6;

    if (two_headers)
        // We need space for the name after the +.
        total_size += name_length;

    // This is the canonical way to create an uninitialized bytestring of given size
    PyObject * retval = PyBytes_FromStringAndSize(NULL, total_size);
    if (retval == NULL)
        return PyErr_NoMemory();

    char * retval_ptr = PyBytes_AS_STRING(retval);

    // Write the sequences into the bytestring at the correct positions.
    size_t cursor;
    retval_ptr[0] = '@';
    memcpy(retval_ptr + 1, name, name_length);
    cursor = name_length + 1;
    retval_ptr[cursor] = '\n'; cursor += 1;
    memcpy(retval_ptr + cursor, sequence, sequence_length);
    cursor += sequence_length;
    retval_ptr[cursor] = '\n'; cursor += 1;
    retval_ptr[cursor] = '+'; cursor += 1;
    if (two_headers){
        memcpy(retval_ptr + cursor, name, name_length);
        cursor += name_length;
    }
    retval_ptr[cursor] = '\n'; cursor += 1;
    memcpy(retval_ptr + cursor, qualities, qualities_length);
    cursor += qualities_length;
    retval_ptr[cursor] = '\n';
    return retval;
}

PyDoc_STRVAR(SequenceRecord_fastq_bytes_two_headers__doc__,
"Return this record in FASTQ format as a bytes object where the header\n"
"(after the @) is repeated on the third line.");

#define SEQUENCE_FASTQ_BYTES_TWO_HEADERS_METHODDEF    \
    {"fastq_bytes_two_headers", \
     (PyCFunction)(void(*)(void))SequenceRecord_fastq_bytes_two_headers, \
     METH_NOARGS, SequenceRecord_fastq_bytes_two_headers__doc__}

static PyObject *
SequenceRecord_fastq_bytes_two_headers(SequenceRecord *self, PyObject *Py_UNUSED(ignore))
{
    Py_INCREF(Py_True);
    PyObject *args[] = {Py_True};
    PyObject * retval = SequenceRecord_fastq_bytes(self, args, (Py_ssize_t)1, NULL);
    Py_DECREF(Py_True);
    return retval;
}

PyDoc_STRVAR(SequenceRecord_qualities_as_bytes__doc__,
"Return the qualities as a bytes object.\n\n"
"This is a faster version of qualities.encode('ascii').");

#define SEQUENCE_QUALITIES_AS_BYTES_METHODDEF    \
    {"qualities_as_bytes", \
     (PyCFunction)(void(*)(void))SequenceRecord_qualities_as_bytes, \
     METH_NOARGS, SequenceRecord_qualities_as_bytes__doc__}

static PyObject *
SequenceRecord_qualities_as_bytes(SequenceRecord *self, PyObject *Py_UNUSED(ignore))
{
    return PyUnicode_AsASCIIString(self->qualities);
}

PyDoc_STRVAR(SequenceRecord_is_mate__doc__,
"Check whether this instance and another are part of the same read pair\n\n");

#define SEQUENCERECORD_IS_MATE_METHODDEF    \
    {"is_mate", \
     (PyCFunction)(void(*)(void))SequenceRecord_is_mate, \
     METH_O, SequenceRecord_is_mate__doc__}

static PyObject *
SequenceRecord_is_mate(SequenceRecord *self, SequenceRecord *other)
{
    if (!(Py_TYPE(other) == &SequenceRecord_Type)) {
        PyErr_SetString(PyExc_TypeError, "other must be a SequenceRecord object.");
        return NULL;
    }
    char * header1_chars = PyUnicode_DATA(self->name);
    size_t header1_length = PyUnicode_GET_LENGTH(self->name);
    char * header2_chars = PyUnicode_DATA(other->name);
    return PyBool_FromLong(
        record_ids_match(header1_chars, header2_chars, header1_length));
}

PyDoc_STRVAR(SequenceRecord_reverse_complement__doc__,
"Return the reverse complement of this record."
);

#define SEQUENCE_REVERSE_COMPLEMENT_METHODDEF \
    {"reverse_complement", (PyCFunction)(void(*)(void))SequenceRecord_reverse_complement, \
     METH_NOARGS, SequenceRecord_reverse_complement__doc__}

static PyObject *
SequenceRecord_reverse_complement(SequenceRecord *self, PyObject *Py_UNUSED(noargs)) 
{
    Py_ssize_t sequence_length = PyUnicode_GET_LENGTH(self->sequence);
    PyObject *reversed_sequence_obj = PyUnicode_New(sequence_length, 127);
    PyObject *reversed_qualities_obj = NULL;
    if (reversed_sequence_obj == NULL) {
        return PyErr_NoMemory();
    }
    char *reversed_sequence = PyUnicode_DATA(reversed_sequence_obj);
    char *sequence = PyUnicode_DATA(self->sequence);
    unsigned char nucleotide;

    Py_ssize_t reverse_cursor = sequence_length;
    Py_ssize_t cursor;
    for (cursor = 0; cursor < sequence_length; cursor += 1) {
        reverse_cursor -= 1;
        nucleotide = sequence[cursor];
        reversed_sequence[reverse_cursor] = NUCLEOTIDE_COMPLEMENTS[nucleotide];
    }
    
    if (self->qualities != NULL) {
        reverse_cursor = sequence_length;
        reversed_qualities_obj = PyUnicode_New(sequence_length, 127);
        if (reversed_qualities_obj == NULL) {
            Py_DECREF(reversed_sequence_obj);
            return PyErr_NoMemory();
        }
        char *reversed_qualities = PyUnicode_DATA(reversed_qualities_obj);
        char *qualities = PyUnicode_DATA(self->qualities);
        for (cursor = 0; cursor < sequence_length; cursor += 1) {
            reverse_cursor -= 1;
            reversed_qualities[reverse_cursor] = qualities[cursor];
        }
    }
    else {
        reversed_qualities_obj = NULL;
    }
    Py_INCREF(self->name);
    return new_sequence_record(self->name, reversed_sequence_obj, reversed_qualities_obj);
    
}

static PyMethodDef SequenceRecord_methods[] = {
    SEQUENCE_FASTQ_BYTES_METHODDEF,
    SEQUENCE_FASTQ_BYTES_TWO_HEADERS_METHODDEF,
    SEQUENCE_QUALITIES_AS_BYTES_METHODDEF,
    SEQUENCERECORD_IS_MATE_METHODDEF,
    SEQUENCE_REVERSE_COMPLEMENT_METHODDEF,
    {NULL}
};


// MAPPING METHODS
static Py_ssize_t 
SequenceRecord__len__(SequenceRecord *self) {
    return PyObject_Size(self->sequence);
}

static PyObject * 
SequenceRecord_get_item(SequenceRecord *self, PyObject *key) 
{
    PyObject * qualities;
    
    PyObject * sequence = PyObject_GetItem(self->sequence, key);
    if (sequence == NULL) {
        return NULL;
    }
    if (self->qualities == NULL) {
        qualities = NULL;
    } else {
        qualities = PyObject_GetItem(self->qualities, key);
        if (qualities == NULL) {
            return NULL;
        }
    }
    Py_INCREF(self->name);
    return new_sequence_record(self->name, sequence, qualities);
}

static PyMappingMethods SequenceRecordMappingMethods = {
    .mp_length = (lenfunc)SequenceRecord__len__,
    .mp_subscript = (binaryfunc)SequenceRecord_get_item,
};

static PyTypeObject SequenceRecord_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_sequence.SequenceRecord",
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_basicsize = sizeof(SequenceRecord),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor)SequenceRecord_dealloc,
    .tp_doc = SequenceRecord__init____doc__,
    .tp_new = SequenceRecord__new__,
    .tp_getset = SequenceRecord_properties,
    .tp_methods = SequenceRecord_methods,
    .tp_repr = (reprfunc)SequenceRecord__repr__,
    .tp_richcompare = (richcmpfunc)SequenceRecord__richcompare__,
    .tp_as_mapping = &SequenceRecordMappingMethods,
};


typedef struct {
    PyObject_HEAD
    Py_ssize_t buffer_size;
    char * buffer;
    Py_ssize_t bytes_in_buffer;
    PyObject * sequence_class;
    int use_custom_class;
    int extra_newline;
    int yielded_two_headers;
    int eof;
    PyObject * file;
    PyObject * read_method;
    char * record_start;
    Py_ssize_t number_of_records;
} FastqIter;

static void 
FastqIter_dealloc(FastqIter *self) {
    Py_CLEAR(self->file);
    Py_CLEAR(self->sequence_class);
    Py_CLEAR(self->read_method);
    PyMem_Free(self->buffer);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyTypeObject FastqIter_Type;

static PyObject *
Fastqiter__new__(PyTypeObject *subtype, PyObject *args, PyObject *kwargs) {
    PyObject *file = NULL;
    PyObject *sequence_class = (PyObject *)&SequenceRecord_Type;
    Py_ssize_t buffer_size = 128 * 1024;

    static char * _keywords[] = {"file", "sequence_class", "buffer_size", NULL};
    static char * _format = "OO!n|:SequenceRecord";
    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, _format, _keywords,
        &file,
        (PyObject *)&PyType_Type, &sequence_class,
        &buffer_size)) {
        return NULL;
    }
    if (buffer_size < 1) {
      PyErr_SetString(PyExc_ValueError, "Starting buffer size too small");
      return NULL;
    }
    FastqIter * self = PyObject_New(FastqIter, subtype);
    self->buffer_size = buffer_size;
    self->buffer = PyMem_Malloc(buffer_size);
    if (self->buffer == NULL) {
      return PyErr_NoMemory();
    }
    self->record_start = self->buffer;
    self->bytes_in_buffer = 0;
    Py_INCREF(sequence_class);
    self->sequence_class = sequence_class;
    self->use_custom_class = (sequence_class != (PyObject *)&SequenceRecord_Type);
    self->number_of_records = 0;
    self->extra_newline = 0; 
    self->yielded_two_headers = 0;
    self->eof = 0;
    Py_INCREF(file);
    self->file = file;
    self->read_method = PyUnicode_FromString("read");
    return (PyObject *)self;
}

static PyObject *FastqFormatError;

static inline void 
raise_FastqFormatError(PyObject *message, Py_ssize_t line) {
    PyObject *line_obj;
    if (line < 0) {
        line_obj = Py_None;
    } else {
        line_obj = PyLong_FromSsize_t(line);
    }
    PyObject *args = PyTuple_New(2);
    PyTuple_SET_ITEM(args, 0, message);
    PyTuple_SET_ITEM(args, 1, line_obj);
    PyObject *err = PyObject_CallObject(FastqFormatError, args);
    Py_DECREF(args);
    PyErr_SetObject(FastqFormatError, err);
}

static int 
FastqIter__read_into_buffer(FastqIter *self) {
    // This function sets self.record_start at 0 and makes sure self.buffer
    // starts at the start of a FASTQ record. Any incomplete FASTQ remainder
    // of the already processed buffer is moved to the start of the buffer
    // and the rest of the buffer is filled up with bytes from the file.
    char * tmp;
    Py_ssize_t remaining_bytes;
    if ((self->record_start == self->buffer) && self->bytes_in_buffer == self->buffer_size) {
      // Buffer too small, double it.
      self->buffer_size *= 2;
      tmp = PyMem_Realloc(self->buffer, self->buffer_size); 
      if (tmp == NULL) {
        PyErr_NoMemory();
        return -1;
      }
      self->buffer = tmp;
    }
    else {
      // Move the incomplete record from the end of the buffer to the beginning.
      remaining_bytes = self->bytes_in_buffer - (self->record_start - self->buffer);
      // Memmove copies safely when dest and src overlap.
      memmove(self->buffer, self->record_start, remaining_bytes);
      self->bytes_in_buffer = remaining_bytes;
    }
    self->record_start = self->buffer;

    Py_ssize_t empty_bytes_in_buffer = self->buffer_size - self->bytes_in_buffer;
    PyObject * filechunk = PyObject_CallMethodObjArgs(
      self->file, self->read_method, PyLong_FromSsize_t(empty_bytes_in_buffer), NULL);
    if (filechunk == NULL || !PyBytes_CheckExact(filechunk)) {
        PyErr_SetString(PyExc_TypeError, "self.file is not a binary file reader.");
        Py_DECREF(filechunk);
        return -1;
    }
    Py_ssize_t filechunk_size = PyBytes_GET_SIZE(filechunk);
    if (filechunk_size > empty_bytes_in_buffer) {
        PyErr_Format(PyExc_ValueError, 
                    "read() returned too much data: %ld bytes requested, "
                    "%ld bytes returned.", empty_bytes_in_buffer, filechunk_size);
        Py_DECREF(filechunk);
        return -1;
    }
    memcpy(self->buffer + self->bytes_in_buffer, 
           PyBytes_AS_STRING(filechunk), filechunk_size);
    Py_DECREF(filechunk);
  
    if (!string_is_ascii(self->buffer + self->bytes_in_buffer, filechunk_size)) {
        raise_FastqFormatError(
            PyUnicode_FromString("Non-ASCII characters found in record."), 
            -1);
        return -1;
    }
    self->bytes_in_buffer += filechunk_size;

    if (filechunk_size == 0) {  // End of file
        if (self->bytes_in_buffer == 0) {
            // All records are processed
            self->eof = 1;
        }
        else if (!self->extra_newline && self->buffer[self->bytes_in_buffer -1] != '\n') {
            // There is still data in the buffer and its last character is
            // not a newline: This is a file that is missing the final 
            // Newline. Append a newline and continue.
            self->buffer[self->bytes_in_buffer] = '\n';
            self->bytes_in_buffer += 1;
            self->extra_newline = 1;
        }
        else { // Incomplete FASTQ records are present.
            if (self->extra_newline) {
                // Do not report the linefeed that was added by dnaio but
                // was not present in the original input.
                self->bytes_in_buffer -= 1;
            }
            PyObject *record = PyUnicode_DecodeASCII(
                        self->record_start, self->bytes_in_buffer, NULL);
            Py_ssize_t record_line_count = PyUnicode_Count(
                record, 
                PyUnicode_FromString("\n"),
                0,
                self->bytes_in_buffer
            );
            raise_FastqFormatError(
                PyUnicode_FromFormat( 
                    "Premature end of file encountered. The incomplete final "
                    "record was: %500S",
                    record
                ),
                self->number_of_records * 4 + record_line_count);
            return -1;
        }
    }
    return 0;
}

static PyObject * 
FastqIter_iter(PyObject * self){
    Py_INCREF(self);
    return self;
}

static PyObject *
FastqIter_next(FastqIter * self) {
    PyObject * retval;
    PyObject * name;
    PyObject * sequence;
    PyObject * qualities;
    char * buffer_end;
    char * name_start;
    char * name_end;
    char * sequence_start; 
    char * sequence_end;
    char * second_header_start;
    char * second_header_end;
    char * qualities_start;
    char * qualities_end; 
    Py_ssize_t name_length, sequence_length, second_header_length, qualities_length;
    // Repeatedly attempt to parse the buffer until we have found a full record.
    // If an attempt fails, we read more data before retrying.
    while (1) {
        buffer_end = self->buffer + self->bytes_in_buffer;
        if (self->eof) {
            PyErr_SetNone(PyExc_StopIteration);
            return NULL;
        }
        name_end = memchr(self->record_start, '\n', (buffer_end - self->record_start));
        if (name_end == NULL) {
            if (FastqIter__read_into_buffer(self) != 0) {
              return NULL;
            }
            continue;
        }
        
        sequence_start = name_end + 1;
        sequence_end = memchr(sequence_start, '\n', (buffer_end - sequence_start));
        if (sequence_end == NULL) {
            if (FastqIter__read_into_buffer(self) != 0) {
                return NULL;
            }
            continue;
        }
        
        second_header_start = sequence_end + 1;
        second_header_end = memchr(second_header_start, '\n', (buffer_end - second_header_start));
        if (second_header_end == NULL) {
            if (FastqIter__read_into_buffer(self) != 0) {
                return NULL;
            }
            continue;
        }

        qualities_start = second_header_end + 1;
        qualities_end = memchr(qualities_start, '\n', (buffer_end - qualities_start));
        if (qualities_end == NULL) {
            if (FastqIter__read_into_buffer(self) != 0) {
                return NULL;
            }
            continue;
        }

        if (self->record_start[0] != '@') {
            raise_FastqFormatError(
                PyUnicode_FromFormat(
                    "Line expected to start with '@' but found '%c'", 
                    self->record_start[0]),
                self->number_of_records * 4
            );
            return NULL;
        }

        if (second_header_start[0] != '+') {
            raise_FastqFormatError(
                PyUnicode_FromFormat( 
                    "Line expected to start with '+' but found '%c'", 
                    second_header_start[0]),
                self->number_of_records * 4 + 2
            );
            return NULL;
        }

        name_start = self->record_start + 1;  // skip @
        second_header_start = second_header_start + 1;  // skip + 
        name_length = name_end - name_start;
        sequence_length = sequence_end - sequence_start;
        second_header_length = second_header_end - second_header_start;
        qualities_length = qualities_end - qualities_start;

        // Check for \r\n line-endings and compensate;
        if ((name_end -1)[0] == '\r') {
            name_length -= 1;
        }
        if ((sequence_end - 1)[0] == '\r') {
            sequence_length -= 1;
        }
        if ((second_header_end - 1)[0] == '\r') {
            second_header_length -= 1;
        }
        if ((qualities_end - 1)[0] == '\r') {
            qualities_length -= 1;
        }

        if (second_header_length) {
            if ((name_length != second_header_length) || memcmp(second_header_start, name_start, second_header_length) !=0) {
                raise_FastqFormatError(
                    PyUnicode_FromFormat(
                        "Sequence descriptions don't match (%R != %R).\n"
                        "The second sequence header must be either empty or equal "
                        "to the first description",
                        PyUnicode_DecodeASCII(name_start, name_length, "strict"),
                        PyUnicode_DecodeASCII(second_header_start, second_header_length, "strict")
                    ),
                    self->number_of_records * 4 + 2
                );
                return NULL;
            }
        }

        if (sequence_length != qualities_length) {
            raise_FastqFormatError(
                PyUnicode_FromString("Length of sequence and qualities differ"),
                self->number_of_records * 4 + 3
            );
            return NULL;
        }

        if ((self->number_of_records == 0) && !(self->yielded_two_headers)) {
            self->yielded_two_headers = 1;
            return PyBool_FromLong(second_header_length);
        }
        
        name = PyUnicode_New(name_length, 127);
        sequence = PyUnicode_New(sequence_length, 127);
        qualities = PyUnicode_New(qualities_length, 127);
        if ((name == NULL) || (sequence == NULL) || (qualities == NULL)) {
            return PyErr_NoMemory();
        }
        memcpy(PyUnicode_DATA(name), name_start, name_length);
        memcpy(PyUnicode_DATA(sequence), sequence_start, sequence_length);
        memcpy(PyUnicode_DATA(qualities), qualities_start, qualities_length);

        if (self->use_custom_class) {
            retval = PyObject_CallFunctionObjArgs(self->sequence_class, name, sequence, qualities);
        }
        else {
            retval = new_sequence_record(name, sequence, qualities);
        }

        self->number_of_records += 1; 
        self->record_start = qualities_end + 1;
        return retval;
    }
}

static PyMemberDef FastqIter_members[] = {
    {"number_of_records", T_PYSSIZET, offsetof(FastqIter, number_of_records), 
     READONLY, NULL},
    {NULL},
};


static PyTypeObject FastqIter_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_core.FastqIter",
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_basicsize = sizeof(FastqIter),
    .tp_dealloc = (destructor)FastqIter_dealloc,
    .tp_new = Fastqiter__new__,
    .tp_iter = FastqIter_iter,
    .tp_iternext = (iternextfunc)FastqIter_next,
    .tp_members = FastqIter_members,
};


PyDoc_STRVAR(paired_fastq_heads__doc__,
"Skip forward in the two buffers by multiples of four lines.\n"
"\n"
"Return a tuple (length1, length2) such that buf1[:length1] and\n"
"buf2[:length2] contain the same number of lines (where the\n"
"line number is divisible by four).\n"
);

#define PAIRED_FASTQ_HEADS_METHODDEF \
    {"paired_fastq_heads", (PyCFunction)(void(*)(void))paired_fastq_heads, \
     METH_VARARGS | METH_KEYWORDS, paired_fastq_heads__doc__}

static PyObject *
paired_fastq_heads(PyObject *module, PyObject *args, PyObject *kwargs)
{
    char *keywords[] = {"buf1", "buf2", "end1", "end2", NULL};
    char *format = "y*y*nn|:paired_fastq_heads";
    Py_buffer buf1;
    Py_buffer buf2;
    Py_ssize_t end1;
    Py_ssize_t end2;
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, format, keywords,
            &buf1, &buf2, &end1, &end2)) {
        return NULL;
    }
    Py_ssize_t linebreaks = 0;
    char *data1 = buf1.buf;
    char *data2 = buf2.buf;
    char *data1_end = data1 + Py_MIN(end1, buf1.len);
    char *data2_end = data2 + Py_MIN(end2, buf2.len);
    char *pos1 = data1;
    char *pos2 = data2;
    char *record_start1 = data1;
    char *record_start2 = data2;

    while (1) {
        pos1 = memchr(pos1, '\n', data1_end - pos1);
        if (pos1 == NULL) 
            break;
        pos1 += 1;
        pos2 = memchr(pos2, '\n', data2_end - pos2);
        if (pos2 == NULL) 
            break;
        pos2 += 1;
        linebreaks += 1;
        if (linebreaks == 4) {
            linebreaks = 0;
            record_start1 = pos1;
            record_start2 = pos2;
        }
    }
    // Hit the end of the data block
    // This code will always be reached, so the buffers are always safely released.
    PyBuffer_Release(&buf1);
    PyBuffer_Release(&buf2);
    PyObject *record_end1 = PyLong_FromSize_t(record_start1 - data1);
    PyObject *record_end2 = PyLong_FromSize_t(record_start2 - data2);
    PyObject *retval = PyTuple_New(2);
    PyTuple_SET_ITEM(retval, 0, record_end1);
    PyTuple_SET_ITEM(retval, 1, record_end2);
    return retval;
}


PyDoc_STRVAR(record_names_match__doc__,
"Check whether the sequence record ids id1 and id2 are compatible, ignoring a\n"
"suffix of '1', '2' or '3'. This exception allows to check some old\n"
"paired-end reads that have IDs ending in '/1' and '/2'. Also, the\n"
"fastq-dump tool (used for converting SRA files to FASTQ) appends '.1', '.2\n"
"and sometimes '.3' to paired-end reads if option -I is used.\n"
"\n"
"\nDeprecated, use `SequenceRecord.is_mate` instead\n");

#define RECORD_NAMES_MATCH_METHODDEF \
    {"record_names_match", (PyCFunction)(void(*)(void))record_names_match, \
     METH_VARARGS | METH_KEYWORDS, record_names_match__doc__}

static PyObject *
record_names_match(PyObject *module, PyObject *args, PyObject *kwargs) 
{
    PyObject *header1 = NULL;
    PyObject *header2 = NULL;
    char *keywords[] = {"header1", "header2", NULL};
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "O!O!|:record_names_match", keywords,
            &PyUnicode_Type, &header1,
            &PyUnicode_Type, &header2)) {
        return NULL;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(header1)) {
        PyErr_SetString(
            PyExc_ValueError, 
            "header1 must be a valid ASCII-string."
        );
        return NULL;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(header2)) {
        PyErr_SetString(
            PyExc_ValueError, 
            "header2 must be a valid ASCII-string."
        );
        return NULL;
    }
    char *header1chars = PyUnicode_DATA(header1);
    char *header2chars = PyUnicode_DATA(header2); 
    size_t header1_length = PyUnicode_GET_LENGTH(header1);
    return PyBool_FromLong(
        record_ids_match(header1chars, header2chars, header1_length));
}


PyDoc_STRVAR(records_are_mates__doc__,
"Check if the provided `SequenceRecord` objects are all mates of each other by\n"
"comparing their record IDs.\n"
"Accepts two or more `SequenceRecord` objects.\n"
"\n"
"This is the same as `SequenceRecord.is_mate` in the case of only two records,\n"
"but allows for for cases where information is split into three records or more\n"
"(such as UMI, R1, R2 or index, R1, R2).\n"
"\n"
"If there are only two records to check, prefer `SequenceRecord.is_mate`.\n"
"\n"
"Example usage::\n"
"\n"
"    for records in zip(*all_my_fastq_readers):\n"
"        if not records_are_mates(*records):\n"
"            raise MateError(f\"IDs do not match for {records}\")\n"
"\n"
"Args:\n"
"    *args: two or more `~dnaio.SequenceRecord` objects\n"
"\n"
"Returns: True or False\n"
);

#define RECORDS_ARE_MATES_METHODDEF \
    {"records_are_mates", (PyCFunction)(void(*)(void))records_are_mates, \
     METH_FASTCALL, records_are_mates__doc__}

static PyObject *
records_are_mates(PyObject *module, PyObject *const *args, Py_ssize_t nargs) {
    if (nargs < 2) {
        PyErr_SetString(
            PyExc_TypeError, 
            "records_are_mates requires at least two arguments");
        return NULL;
    }
    PyTypeObject *seqrecord_type = (PyTypeObject *)(&SequenceRecord_Type);

    PyObject *first = args[0];
    if (Py_TYPE(first) != seqrecord_type) {
        PyErr_Format(
            PyExc_TypeError, 
            "%R is not a SequenceRecord object",
            first);
        return NULL;
    }
    SequenceRecord *first_seq = (SequenceRecord *)first;
    PyObject *first_name_obj = first_seq->name;
    char *first_name = PyUnicode_DATA(first_name_obj);
    size_t first_id_length = strcspn(first_name, " \t");
    char first_id_end = first_name[first_id_length - 1];
    int first_id_ends_with_number = first_id_end >= '1' && first_id_end <= '3';
    int are_mates = 1;

    for (Py_ssize_t i=1; i<nargs; i++) {
        PyObject *other = args[i];
        if (Py_TYPE(other) != seqrecord_type) {
            PyErr_Format(
                PyExc_TypeError, 
                "%R is not a SequenceRecord object",
                other);
        return NULL;
        }
        SequenceRecord *other_seq = (SequenceRecord *)other;
        PyObject *other_name_obj = other_seq->name;
        char *other_name = PyUnicode_DATA(other_name_obj);
        size_t other_name_length = PyUnicode_GET_LENGTH(other_name_obj);
        are_mates &= record_ids_match_partial(
            first_name, other_name, first_id_length, other_name_length, 
            first_id_ends_with_number);
    }
    return PyBool_FromLong(are_mates);
}


static PyMethodDef _core_methods[] = {
    RECORD_NAMES_MATCH_METHODDEF,
    RECORDS_ARE_MATES_METHODDEF,
    PAIRED_FASTQ_HEADS_METHODDEF,
    {NULL},
};

static struct PyModuleDef _core_module = {
    PyModuleDef_HEAD_INIT,
    "_core",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,
    _core_methods  /* module methods */
};



PyMODINIT_FUNC
PyInit__core(void)
{
    PyObject *m;

    PyObject *exceptions_module = PyImport_ImportModule("dnaio.exceptions");
    if (exceptions_module == NULL) {
        return NULL;
    }
    FastqFormatError = PyObject_GetAttrString(exceptions_module, "FastqFormatError");

    m = PyModule_Create(&_core_module);
    if (m == NULL)
        return NULL;
    PyTypeObject * SequenceRecordType = &SequenceRecord_Type;
    if (PyType_Ready(SequenceRecordType) != 0) { 
        return NULL;
    }
    Py_INCREF((PyObject *)SequenceRecordType);
    if (PyModule_AddObject(
            m, "SequenceRecord", (PyObject *)SequenceRecordType) != 0) {
        return NULL;
    }
    // Add aliases for backwards compatibility
    Py_INCREF((PyObject *)SequenceRecordType);
    if (PyModule_AddObject(
            m, "Sequence", (PyObject *)SequenceRecordType) != 0) {
        return NULL;
    }

    PyTypeObject * FastqIterType = &FastqIter_Type;
    if (PyType_Ready(FastqIterType) != 0) {
        return NULL;
    }
    Py_INCREF((PyObject *)FastqIterType);
    if (PyModule_AddObject(m, "FastqIter", (PyObject *)FastqIterType) !=0) {
        return NULL;
    }

    return m;
}
