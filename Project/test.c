#include <Python.h>
#include <stdlib.h>
#include <math.h>

static PyObject *
newton(PyObject *self, PyObject *args)
{
    PyObject *cb1, *cb2;
	double eps, xn, x1, x0;
	double res1, res2;

    if (!PyArg_ParseTuple(args, "OOdd", &cb1, &cb2, &xn, &eps))
        return 0;

    if (!PyCallable_Check(cb1) || !PyCallable_Check(cb2)) {
        PyErr_SetString(PyExc_TypeError, "newton: a callable is required");
        return 0;
    }

	PyObject *arg = Py_BuildValue("(d)", xn);
	res1 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg));
	res2 = PyFloat_AsDouble(PyObject_CallObject(cb2, arg));
	x1 = xn - res1/res2;
	x0 = xn;
	while(abs(x0-x1)>eps) {
		x0 = x1;
		arg = Py_BuildValue("(d)", x1);
		res1 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg));
		res2 = PyFloat_AsDouble(PyObject_CallObject(cb2, arg));
		x1 = x1 - res1/res2;
	}
	return Py_BuildValue("d", x1);
}

static PyObject *
chebishev(PyObject *self, PyObject *args)
{
    PyObject *cb1, *cb2, *cb3;
	double eps, xn, x1, x0;
	double res1, res2, res3;

    if (!PyArg_ParseTuple(args, "OOOdd", &cb1, &cb2, &cb3, &xn, &eps))
        return 0;

    if (!PyCallable_Check(cb1) || !PyCallable_Check(cb2) || !PyCallable_Check(cb3)) {
        PyErr_SetString(PyExc_TypeError, "chebishev: a callable is required");
        return 0;
    }

	PyObject *arg = Py_BuildValue("(d)", xn);
	res1 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg));
	res2 = PyFloat_AsDouble(PyObject_CallObject(cb2, arg));
	res3 = PyFloat_AsDouble(PyObject_CallObject(cb3, arg));
	x1  = xn - res1/res2-(pow(res1,2)*res3)/(2*pow(res2,3));
	x0 = xn;
	while(abs(x0-x1)>eps) {
		x0 = x1;
		arg = Py_BuildValue("(d)", x1);
		res1 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg));
		res2 = PyFloat_AsDouble(PyObject_CallObject(cb2, arg));
		res3 = PyFloat_AsDouble(PyObject_CallObject(cb3, arg));
		x1  = x1 - res1/res2-(pow(res1,2)*res3)/(2*pow(res2,3));
	}
	return Py_BuildValue("d", x1);
}

static PyObject *
helley(PyObject *self, PyObject *args)
{
    PyObject *cb1, *cb2, *cb3;
	double eps, xn, x1, x0;
	double res1, res2, res3;

    if (!PyArg_ParseTuple(args, "OOOdd", &cb1, &cb2, &cb3, &xn, &eps))
        return 0;

    if (!PyCallable_Check(cb1) || !PyCallable_Check(cb2) || !PyCallable_Check(cb3)) {
        PyErr_SetString(PyExc_TypeError, "helley: a callable is required");
        return 0;
    }

	PyObject *arg = Py_BuildValue("(d)", xn);
	res1 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg));
	res2 = PyFloat_AsDouble(PyObject_CallObject(cb2, arg));
	res3 = PyFloat_AsDouble(PyObject_CallObject(cb3, arg));
	x1  = xn - (2*res1*res2/(2*pow(res2,2)-res1*res3));
	x0 = xn;
	while(abs(x0-x1)>eps) {
		x0 = x1;
		arg = Py_BuildValue("(d)", x1);
		res1 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg));
		res2 = PyFloat_AsDouble(PyObject_CallObject(cb2, arg));
		res3 = PyFloat_AsDouble(PyObject_CallObject(cb3, arg));
		x1  = x1 - (2*res1*res2/(2*pow(res2,2)-res1*res3));
	}
	return Py_BuildValue("d", x1);
}

static PyObject *
inverseInterpolation(PyObject *self, PyObject *args)
{
    PyObject *cb1;
	double eps, xn, xnm1, xnm2, x1, x0, xm1, xm2;
	double res1, res2, res3;

    if (!PyArg_ParseTuple(args, "Odddd", &cb1, &xn, &xnm1, &xnm2, &eps))
        return 0;

    if (!PyCallable_Check(cb1)) {
        PyErr_SetString(PyExc_TypeError, "inverseInterpolation: a callable is required");
        return 0;
    }

	PyObject *arg = Py_BuildValue("(d)", xn);
	res1 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg));
	PyObject *arg1 = Py_BuildValue("(d)", xnm1);
	res2 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg1));
	PyObject *arg2 = Py_BuildValue("(d)", xnm2);
	res3 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg2));
	x1  = res2*res1*xnm2/((res3-res2)*(res2-res1)) + res3*res1*xnm1/((res2-res3)*(res2-res1))+res3*res2*xn/((res1-res3)*(res1-res2));
	x0 = xn;
	xm1 = xnm1;
	xm2 = xnm2;
	while(abs(x0-x1)>eps) {
        xm2 = xm1;
        arg = Py_BuildValue("(d)", xm1);
		res3 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg));
		xm1 = x0;
		arg1 = Py_BuildValue("(d)", x0);
		res2 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg1));
		x0 = x1;
		arg2 = Py_BuildValue("(d)", x1);
		res1 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg2));
        x1 = res2*res1*xm2/((res3-res2)*(res3-res1)) + res3*res1*xm1/((res2-res3)*(res2-res1))+res3*res2*x0/((res1-res3)*(res1-res2));
	}
	return Py_BuildValue("d", x1);
}

static PyObject *
secant(PyObject *self, PyObject *args)
{
    PyObject *cb1;
	double eps, x_prev, x_curr, tmp, x_next=0;
	double res1, res2;

    if (!PyArg_ParseTuple(args, "Oddd", &cb1, &x_prev, &x_curr, &eps))
        return 0;

    if (!PyCallable_Check(cb1)) {
        PyErr_SetString(PyExc_TypeError, "secant: a callable is required");
        return 0;
    }

	PyObject *arg = Py_BuildValue("(d)", x_prev);
	res1 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg));
	PyObject *arg1 = Py_BuildValue("(d)", x_curr);
	res2 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg1));
	while(abs(x_next - x_curr) > eps)
	{
        tmp = x_next;
        x_next = x_curr - res2 * (x_prev - x_curr) / (res1- res2);
        x_prev = x_curr;
        PyObject *arg = Py_BuildValue("(d)", x_prev);
        res1 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg));
        x_curr = tmp;
        PyObject *arg1 = Py_BuildValue("(d)", x_curr);
        res2 = PyFloat_AsDouble(PyObject_CallObject(cb1, arg1));
	}
	return Py_BuildValue("d", x_next);
}

static PyMethodDef ownmod_methods[] = {
    {
        "newton",
        newton,
        METH_VARARGS,
        "newton function"
    },
    {
        "chebishev",
        chebishev,
        METH_VARARGS,
        "chebishev function"
    },
    {
        "helley",
        helley,
        METH_VARARGS,
        "helley function"
    },
    {
        "inverseInterpolation",
        inverseInterpolation,
        METH_VARARGS,
        "inverseInterpolation function"
    },
    {
        "secant",
        secant,
        METH_VARARGS,
        "secant function"
    },
    { NULL, NULL, 0, NULL }
};


static PyModuleDef simple_module = {
    /* Info about module */
    PyModuleDef_HEAD_INIT,
    "methods", // my_plus.__name__
    "amazing documentation", // my_plus.__doc__
    -1,
    ownmod_methods, // methods are here
    NULL,
    NULL,
    NULL,
    NULL
};


PyMODINIT_FUNC PyInit_methods(void)
{
    PyObject* m;
    // creating module object
    m = PyModule_Create(&simple_module);
    if (m == NULL)
        return NULL;

    return m;
}
