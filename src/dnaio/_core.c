#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"         // PyMemberDef

#define PY_SSIZE_T_CLEAN
#include <Python.h>

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

static inline PyObject *
create_fastq_record(char * name, char * sequence, char * qualities,
                    Py_ssize_t name_length,
                    Py_ssize_t sequence_length,
                    Py_ssize_t qualities_length,
                    int two_headers) {
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
        &qualities))
        return NULL;
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

static inline PyObject * 
sequence_to_fastq_record_impl(SequenceRecord *self, int two_headers){
    if (self->qualities == NULL) {
        PyErr_SetString(PyExc_ValueError, 
        "Cannot create FASTQ bytes from a sequence without qualities.");
    }
    Py_ssize_t name_length = PyUnicode_GetLength(self->name);
    Py_ssize_t sequence_length = PyUnicode_GetLength(self->sequence);
    Py_ssize_t qualities_length = PyUnicode_GetLength(self->qualities);
    if (name_length == -1 || sequence_length == -1 || qualities_length == -1)
            // Not of type PyUnicode
            return NULL;
    if (!(
        (PyUnicode_KIND(self->name) == PyUnicode_1BYTE_KIND) && 
        (PyUnicode_KIND(self->sequence)== PyUnicode_1BYTE_KIND) &&
        (PyUnicode_KIND(self->qualities) == PyUnicode_1BYTE_KIND))) {
            PyErr_SetString(
                PyExc_ValueError, 
                "Name, sequence and qualities must all be valid ASCII strings."
            );
    }
    // Unsafe macros as type is already checked.
    char * name = (char *)PyUnicode_1BYTE_DATA(self->name);
    char * sequence = (char *)PyUnicode_1BYTE_DATA(self->sequence);
    char * qualities = (char *)PyUnicode_1BYTE_DATA(self->qualities);
    return create_fastq_record(name, sequence, qualities,
                               name_length, sequence_length, qualities_length,
                               two_headers);
}

PyDoc_STRVAR(SequenceRecord_fastq_bytes__doc__,
"Return the entire FASTQ record as bytes which can be written\n"
"into a file.");

#define SEQUENCE_FASTQ_BYTES_METHODDEF    \
    {"fastq_bytes", (PyCFunction)(void(*)(void))SequenceRecord_fastq_bytes, \
     METH_NOARGS, SequenceRecord_fastq_bytes__doc__}

static PyObject *
SequenceRecord_fastq_bytes(SequenceRecord *self, PyObject *Py_UNUSED(ignore)){
    return sequence_to_fastq_record_impl(self, 0);
}

#define BYTES_SEQUENCE_FASTQ_BYTES_METHODDEF    \
    {"fastq_bytes", (PyCFunction)(void(*)(void))BytesSequenceRecord_fastq_bytes, \
    METH_NOARGS, SequenceRecord_fastq_bytes__doc__}

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
    return sequence_to_fastq_record_impl(self, 1);
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

static PyMethodDef SequenceRecord_methods[] = {
    SEQUENCE_FASTQ_BYTES_METHODDEF,
    SEQUENCE_FASTQ_BYTES_TWO_HEADERS_METHODDEF,
    SEQUENCE_QUALITIES_AS_BYTES_METHODDEF,
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

static PyTypeObject SequenceRecord_type = {
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


static struct PyModuleDef _sequence_module = {
    PyModuleDef_HEAD_INIT,
    "_sequence",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,
    NULL  /* module methods */
};


PyMODINIT_FUNC
PyInit__sequence(void)
{
    PyObject *m;

    m = PyModule_Create(&_sequence_module);
    if (m == NULL)
        return NULL;
    PyTypeObject * SequenceRecordType = &SequenceRecord_type;
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
    return m;
}
