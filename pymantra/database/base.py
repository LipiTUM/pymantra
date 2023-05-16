import traceback
from typing import Union, Tuple
from neo4j import GraphDatabase, Session, Driver, Result
from neo4j.exceptions import ServiceUnavailable


class Neo4jBaseConnector:
    """
    Base class to connect to a neo4j database with utilities to check
    connectivity and use a class in a with... statement

    Attributes
    ----------
    driver : neo4j.Driver
        database driver

    session : neo4j.Session
        database session
    """
    driver: Driver
    session: Session

    def __init__(
        self, uri: str, auth: Union[Tuple[str, str], None] = None, **kwargs
    ):
        """
        Parameters
        ----------
        uri : str
            database uri

        auth : Tuple[str, str], Optional
            Optional database credentials int the form of (user, password).
            If database is not secured pass None.

        data : ResourceData (NamedTuple)
            Data to dump into the database

        verbose : bool, Optional, default True
            Whether function should be verbose

        reset_at_start : bool, Optional, default False
            If True the database will be emptied before data is added

        Raises
        ------
        ConnectionError
            If database is not reachable with the given parameters
        """
        self.driver = GraphDatabase.driver(uri, auth=auth, **kwargs)
        self.session = self.driver.session()
        if not self.active_connection():
            if auth is None:
                add = "If your database requires login credentials please " \
                      "add them via the 'auth' parameter"
            else:
                add = "Please make sure you parsed the correct login " \
                      "credentials"
            raise ConnectionError(
                f"Connection to {uri} failed. {add}. If you have not "
                "installed mantra's neo4j database locally or if you are "
                "running the local database in the pymantra-db-api docker "
                "please use the 'APINetworkGenerator' class for all graph "
                "queries."
            )

    def valid_port(self) -> bool:
        try:
            # this is just a simple query to check the driver connection
            self.session.run(
                "MATCH (n:metabolite) RETURN count(n) as count")
            return True
        except ServiceUnavailable:
            return False

    def active_connection(self) -> bool:
        """
        Verifying the connection to the database

        Returns
        -------
        bool
            True if connection is established, else False
        """
        try:
            # this is just a simple query to check the driver connection
            self.session.run(
                "MATCH (n:metabolite) RETURN count(n) as count")
            return True
        except ServiceUnavailable:
            return False

    def close(self):
        """
        Closing the database connection. Should always be called
        once the object instance is not needed anymore
        """
        if hasattr(self, "session"):
            # if constructor fails at driver call session is not set
            self.session.close()
        if hasattr(self, "driver"):
            self.driver.close()

    def __del__(self):
        self.close()

    # for usage with with(NetworkGenerator()) as ...
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        self.close()
        if exc_type is not None:
            traceback.print_exception(exc_type, exc_value, tb)

    def run(self, query: str, *parameters, **kwargs) -> Result:
        """Run a neo4j query

        Parameters
        ----------
        query: str
            Query string
        parameters: Dict[str, any], optional
           Query parameters to insert into query
        kwargs
           Query parameters to insert into query, Preceding over those in
           `parameters`

        Returns
        -------
        neo4j.Result
            Query result
        """
        return self.session.run(query, *parameters, **kwargs)
